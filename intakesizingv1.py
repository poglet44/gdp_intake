#!/usr/bin/env python3
"""
prelim_intake_sizing_v7.py

Preliminary supersonic intake sizing tool with rigid station handling
and explicit GasTurb area interpretation.

PURPOSE
-------
This script performs first-pass intake sizing for a supersonic aircraft using:
- Station 0 freestream static inputs
- Capture sizing from actual mass flow
- A rigid station-2 solve
- Explicit interpretation of the user-supplied GasTurb area

KEY NEW FEATURE
---------------
You can now enter:
- a raw station-2 area
- a hub/tip ratio Rh/Rt
- an area interpretation:
    * "gross"   -> input area is full circular face area = pi * Rt^2
    * "annulus" -> input area is annulus flow area = pi * (Rt^2 - Rh^2)

The script will then:
- solve Rt
- solve Rh
- solve D2
- solve station-2 Mach from mdot_2, Tt2, Pt2 and the chosen flow area
- compare all of these quantities clearly in the output

IMPORTANT
---------
This is still a low-order preliminary design tool.
It does NOT yet predict:
- shock system
- throat area
- diffuser losses
- bleed / bypass
- distortion
- off-design behaviour
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional


# =============================================================================
# ASSUMPTIONS USED IN THIS VERSION
# =============================================================================
# 1. Single design-point preliminary sizing only.
#
# 2. Capture area is computed from actual mass flow using station-0 freestream:
#       A_capture = mdot / (rho0 * V0)
#
# 3. Station-0 total conditions are derived from station-0 static conditions
#    using perfect-gas compressible-flow relations.
#
# 4. Station-2 is solved using ONE rigid solve mode at a time.
#
# 5. The user-supplied station-2 area is explicitly interpreted as either:
#       - gross face area
#       - annulus area
#    The script does not mix these silently.
#
# 6. Rh/Rt is treated as an input geometric ratio.
#
# 7. Perfect-gas behaviour is assumed with constant gamma and R.
#
# 8. cp values are not used in this version.
#
# 9. No shock, throat, diffuser-loss, bleed, bypass, distortion, or off-design
#    modelling is included in this version.
#
# 10. 2D bifurcated and axisymmetric concepts are assumed to require
#     approximately the same capture area for the same mdot and freestream state.
# =============================================================================


@dataclass
class IntakeSizingInputs:
    # -------------------------
    # Station 0: freestream static inputs
    # -------------------------
    mach_0: float = 1.8
    altitude_ft: float = 60000.0
    t_0_k: float = 216.65
    p_0_kpa: float = 7.172

    # Perfect gas properties
    gamma: float = 1.4
    r_j_per_kgk: float = 287.05

    # -------------------------
    # Engine / intake inputs
    # -------------------------
    engines_total: int = 4
    mdot_per_engine_kg_s: float = 159.19

    # Practical capture sizing margin
    capture_area_margin_fraction: float = 0.08

    # Example 2D capture widths
    example_2d_widths_m: List[float] = field(
        default_factory=lambda: [2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
    )

    # -------------------------
    # Station 2 solve mode
    # Allowed values:
    #   "solve_from_mdot_Tt2_Pt2_A2"
    #   "solve_from_mdot_Tt2_Pt2_M2"
    # -------------------------
    station2_mode: str = "solve_from_mdot_Tt2_Pt2_A2"

    # -------------------------
    # Station 2 independent thermodynamic inputs
    # -------------------------
    mdot_2_kg_s: float = 159.19
    tt_2_k: float = 357.049
    pt_2_kpa: float = 39.169

    # -------------------------
    # Station 2 area / geometry interpretation
    # area_type options:
    #   "gross"
    #   "annulus"
    # -------------------------
    area_2_input_m2: float = 2.020917
    area_2_input_type: str = "gross"
    rh_over_rt: float = 0.45

    # For mode B only
    mach_2_target: float = 0.748

    # -------------------------
    # Optional GasTurb reference values for checking only
    # -------------------------
    ref_mach_2: Optional[float] = 0.748
    ref_t_2_k: Optional[float] = 321.264
    ref_p_2_kpa: Optional[float] = 27.0315
    ref_v_2_m_s: Optional[float] = 268.661
    ref_rho_2_kg_m3: Optional[float] = 0.293123

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Single design-point preliminary sizing only.",
            "Capture area is computed from actual mass flow using station-0 freestream static conditions.",
            "Station-0 total conditions are derived from station-0 static conditions.",
            "Station-2 is solved using exactly one rigid solve mode at a time.",
            "The user-supplied station-2 area is explicitly interpreted as either gross or annulus area.",
            "Rh/Rt is treated as a direct geometry input.",
            "Perfect-gas model with constant gamma and gas constant R.",
            "cp values are not used in this version.",
            "No shock, throat, diffuser-loss, bleed, bypass, distortion, or off-design modelling is included.",
            "2D bifurcated and axisymmetric concepts are assumed to need approximately the same capture area for the same mdot and freestream state.",
        ]
    )


@dataclass
class IntakeSizingResults:
    # Station 0 static
    rho_0_kg_m3: Optional[float]
    a_0_m_s: Optional[float]
    v_0_m_s: Optional[float]
    mass_flux_0_kg_m2_s: Optional[float]

    # Station 0 total
    tt_0_k: Optional[float]
    pt_0_kpa: Optional[float]

    # Capture sizing
    a_capture_ideal_m2: Optional[float]
    a_capture_working_m2: Optional[float]
    d_capture_axisym_ideal_m: Optional[float]
    d_capture_axisym_working_m: Optional[float]

    # Station 2 area interpretation
    area_2_input_m2: Optional[float]
    area_2_input_type: str
    rh_over_rt: Optional[float]
    rt_2_m: Optional[float]
    rh_2_m: Optional[float]
    d2_outer_m: Optional[float]
    a2_gross_m2: Optional[float]
    a2_annulus_m2: Optional[float]
    a2_flow_used_m2: Optional[float]

    # Station 2 solved state
    station2_mode: str
    mdot_2_kg_s: Optional[float]
    tt_2_k: Optional[float]
    pt_2_kpa: Optional[float]
    mach_2: Optional[float]
    t_2_k: Optional[float]
    p_2_kpa: Optional[float]
    rho_2_kg_m3: Optional[float]
    a_2_m_s: Optional[float]
    v_2_m_s: Optional[float]

    # Ratios
    capture_to_a2flow_ideal: Optional[float]
    capture_to_a2flow_working: Optional[float]

    # 2D realised openings
    two_d_examples: List[tuple[float, Optional[float], Optional[float]]]

    # Reference check residuals
    delta_mach_2: Optional[float]
    delta_t_2_k: Optional[float]
    delta_p_2_kpa: Optional[float]
    delta_v_2_m_s: Optional[float]
    delta_rho_2_kg_m3: Optional[float]

    # Warnings
    warnings: List[str]


def circle_area(diameter_m: float) -> float:
    return math.pi * diameter_m**2 / 4.0


def equivalent_diameter_from_area(area_m2: float) -> float:
    return math.sqrt(4.0 * area_m2 / math.pi)


def compute_density(p_kpa: float, t_k: float, r: float) -> float:
    return (p_kpa * 1000.0) / (r * t_k)


def compute_speed_of_sound(gamma: float, r: float, t_k: float) -> float:
    return math.sqrt(gamma * r * t_k)


def total_temperature_from_static(t_static: float, mach: float, gamma: float) -> float:
    return t_static * (1.0 + 0.5 * (gamma - 1.0) * mach**2)


def total_pressure_from_static(p_static: float, mach: float, gamma: float) -> float:
    return p_static * (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (gamma / (gamma - 1.0))


def static_temperature_from_total(tt: float, mach: float, gamma: float) -> float:
    return tt / (1.0 + 0.5 * (gamma - 1.0) * mach**2)


def static_pressure_from_total(pt: float, mach: float, gamma: float) -> float:
    return pt / (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (gamma / (gamma - 1.0))


def mass_flow_parameter(mach: float, gamma: float) -> float:
    return mach * (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (
        -(gamma + 1.0) / (2.0 * (gamma - 1.0))
    )


def mdot_from_total_conditions(
    mach: float,
    area_m2: float,
    pt_kpa: float,
    tt_k: float,
    gamma: float,
    r: float,
) -> float:
    p_pa = pt_kpa * 1000.0
    return (
        area_m2
        * p_pa
        / math.sqrt(tt_k)
        * math.sqrt(gamma / r)
        * mass_flow_parameter(mach, gamma)
    )


def area_required_for_mdot(
    mdot_kg_s: float,
    mach: float,
    pt_kpa: float,
    tt_k: float,
    gamma: float,
    r: float,
) -> float:
    denom = (
        pt_kpa * 1000.0
        / math.sqrt(tt_k)
        * math.sqrt(gamma / r)
        * mass_flow_parameter(mach, gamma)
    )
    return mdot_kg_s / denom


def solve_mach_for_mdot(
    target_mdot: float,
    area_m2: float,
    pt_kpa: float,
    tt_k: float,
    gamma: float,
    r: float,
) -> tuple[Optional[float], Optional[str]]:
    """
    Solve subsonic Mach number for fixed mdot, area, Pt, and Tt.
    Returns (mach, warning_message).
    """
    def f(m: float) -> float:
        return mdot_from_total_conditions(m, area_m2, pt_kpa, tt_k, gamma, r) - target_mdot

    mach_low = 1e-6
    mach_high = 0.999
    tol = 1e-10
    max_iter = 200

    f_low = f(mach_low)
    f_high = f(mach_high)

    if f_low * f_high > 0:
        max_subsonic_mdot = mdot_from_total_conditions(
            mach_high, area_m2, pt_kpa, tt_k, gamma, r
        )
        return None, (
            "Could not bracket a subsonic Mach solution for station 2. "
            f"Demanded mdot = {target_mdot:.4f} kg/s, "
            f"max subsonic mdot under current assumptions ≈ {max_subsonic_mdot:.4f} kg/s."
        )

    low = mach_low
    high = mach_high

    for _ in range(max_iter):
        mid = 0.5 * (low + high)
        f_mid = f(mid)

        if abs(f_mid) < tol:
            return mid, None

        if f_low * f_mid < 0:
            high = mid
        else:
            low = mid
            f_low = f_mid

    return 0.5 * (low + high), "Subsonic Mach solve reached iteration limit; returning last estimate."


def fmt(value: Optional[float], ndp: int = 4) -> str:
    if value is None:
        return "N/A"
    return f"{value:.{ndp}f}"


def derive_geometry_from_area_and_ratio(
    area_input_m2: Optional[float],
    area_type: str,
    rh_over_rt: Optional[float],
    warnings: List[str],
) -> tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Returns:
        rt, rh, d2, a2_gross, a2_annulus, a2_flow_used

    IMPORTANT:
    - a2_gross   = full circular face area
    - a2_annulus = actual flow-passing annulus area
    - a2_flow_used = area used in the station-2 thermodynamic solve

    In this script, the flow solve should always use annulus area.
    """
    rt = rh = d2 = a2_gross = a2_annulus = a2_flow_used = None

    if area_input_m2 is None or area_input_m2 <= 0:
        warnings.append("Station-2 area input must be > 0.")
        return rt, rh, d2, a2_gross, a2_annulus, a2_flow_used

    if rh_over_rt is None or rh_over_rt < 0 or rh_over_rt >= 1:
        warnings.append("Rh/Rt must be within [0, 1).")
        return rt, rh, d2, a2_gross, a2_annulus, a2_flow_used

    area_type_clean = area_type.strip().lower()

    if area_type_clean == "gross":
        # User input is full circular face area
        a2_gross = area_input_m2
        rt = math.sqrt(a2_gross / math.pi)
        rh = rh_over_rt * rt
        d2 = 2.0 * rt
        a2_annulus = math.pi * (rt**2 - rh**2)

        # CRITICAL FIX:
        # even though the input was gross area, the flow passes through the annulus
        a2_flow_used = a2_annulus

    elif area_type_clean == "annulus":
        # User input is annulus flow area
        a2_annulus = area_input_m2
        rt = math.sqrt(a2_annulus / (math.pi * (1.0 - rh_over_rt**2)))
        rh = rh_over_rt * rt
        d2 = 2.0 * rt
        a2_gross = math.pi * rt**2

        # Flow still uses annulus area
        a2_flow_used = a2_annulus

    else:
        warnings.append("area_2_input_type must be either 'gross' or 'annulus'.")
        return rt, rh, d2, a2_gross, a2_annulus, a2_flow_used

    return rt, rh, d2, a2_gross, a2_annulus, a2_flow_used


def compute_intake_sizing(inputs: IntakeSizingInputs) -> IntakeSizingResults:
    warnings: List[str] = []

    # -------------------------------------------------------------------------
    # Basic input checks
    # -------------------------------------------------------------------------
    if inputs.t_0_k <= 0:
        warnings.append("Station-0 static temperature must be > 0 K.")
    if inputs.p_0_kpa <= 0:
        warnings.append("Station-0 static pressure must be > 0 kPa.")
    if inputs.mdot_per_engine_kg_s <= 0:
        warnings.append("Per-engine actual mass flow for capture sizing must be > 0 kg/s.")
    if inputs.tt_2_k <= 0:
        warnings.append("Station-2 total temperature must be > 0 K.")
    if inputs.pt_2_kpa <= 0:
        warnings.append("Station-2 total pressure must be > 0 kPa.")
    if inputs.mdot_2_kg_s <= 0:
        warnings.append("Station-2 mass flow must be > 0 kg/s.")
    if inputs.capture_area_margin_fraction < 0:
        warnings.append("Capture area margin is negative. This may be intentional, but is unusual.")

    # -------------------------------------------------------------------------
    # Station 0 static
    # -------------------------------------------------------------------------
    rho_0 = a_0 = v_0 = mass_flux_0 = None
    if inputs.t_0_k > 0 and inputs.p_0_kpa > 0:
        try:
            rho_0 = compute_density(inputs.p_0_kpa, inputs.t_0_k, inputs.r_j_per_kgk)
            a_0 = compute_speed_of_sound(inputs.gamma, inputs.r_j_per_kgk, inputs.t_0_k)
            v_0 = inputs.mach_0 * a_0
            mass_flux_0 = rho_0 * v_0
        except Exception as e:
            warnings.append(f"Failed to compute station-0 static properties: {e}")

    # -------------------------------------------------------------------------
    # Station 0 total
    # -------------------------------------------------------------------------
    tt_0 = pt_0 = None
    if inputs.t_0_k > 0 and inputs.p_0_kpa > 0:
        try:
            tt_0 = total_temperature_from_static(inputs.t_0_k, inputs.mach_0, inputs.gamma)
            pt_0 = total_pressure_from_static(inputs.p_0_kpa, inputs.mach_0, inputs.gamma)
        except Exception as e:
            warnings.append(f"Failed to compute station-0 total properties: {e}")

    # -------------------------------------------------------------------------
    # Capture sizing
    # -------------------------------------------------------------------------
    a_capture_ideal = a_capture_working = None
    d_capture_axisym_ideal = d_capture_axisym_working = None
    two_d_examples: List[tuple[float, Optional[float], Optional[float]]] = []

    if mass_flux_0 is not None and mass_flux_0 > 0 and inputs.mdot_per_engine_kg_s > 0:
        try:
            a_capture_ideal = inputs.mdot_per_engine_kg_s / mass_flux_0
            a_capture_working = a_capture_ideal * (1.0 + inputs.capture_area_margin_fraction)

            if a_capture_ideal > 0:
                d_capture_axisym_ideal = equivalent_diameter_from_area(a_capture_ideal)
            else:
                warnings.append("Ideal capture area is non-positive.")

            if a_capture_working is not None and a_capture_working > 0:
                d_capture_axisym_working = equivalent_diameter_from_area(a_capture_working)
            else:
                warnings.append("Working capture area is non-positive.")

            for width in inputs.example_2d_widths_m:
                if width > 0:
                    ideal_height = a_capture_ideal / width if a_capture_ideal is not None else None
                    working_height = a_capture_working / width if a_capture_working is not None else None
                else:
                    warnings.append(f"Skipped non-positive 2D width value: {width}")
                    ideal_height = None
                    working_height = None
                two_d_examples.append((width, ideal_height, working_height))

        except Exception as e:
            warnings.append(f"Failed to compute capture sizing: {e}")
    else:
        warnings.append("Capture sizing not performed because station-0 mass flux was invalid.")

    # -------------------------------------------------------------------------
    # Station 2 geometry / area interpretation
    # -------------------------------------------------------------------------
    rt_2 = rh_2 = d2_outer = a2_gross = a2_annulus = a2_flow_used = None
    rt_2, rh_2, d2_outer, a2_gross, a2_annulus, a2_flow_used = derive_geometry_from_area_and_ratio(
        area_input_m2=inputs.area_2_input_m2,
        area_type=inputs.area_2_input_type,
        rh_over_rt=inputs.rh_over_rt,
        warnings=warnings,
    )

    # -------------------------------------------------------------------------
    # Station 2 solve
    # -------------------------------------------------------------------------
    mode = inputs.station2_mode.strip()
    mdot_2 = inputs.mdot_2_kg_s if inputs.mdot_2_kg_s > 0 else None
    tt_2 = inputs.tt_2_k if inputs.tt_2_k > 0 else None
    pt_2 = inputs.pt_2_kpa if inputs.pt_2_kpa > 0 else None

    mach_2 = t_2 = p_2 = rho_2 = a_2_sound = v_2 = None

    if mode == "solve_from_mdot_Tt2_Pt2_A2":
        if a2_flow_used is None or a2_flow_used <= 0:
            warnings.append("Station-2 mode A selected, but the interpreted flow area is invalid.")
        elif mdot_2 is None or tt_2 is None or pt_2 is None:
            warnings.append("Station-2 mode A could not run because mdot_2, Tt2, or Pt2 was invalid.")
        else:
            mach_2, warning = solve_mach_for_mdot(
                target_mdot=mdot_2,
                area_m2=a2_flow_used,
                pt_kpa=pt_2,
                tt_k=tt_2,
                gamma=inputs.gamma,
                r=inputs.r_j_per_kgk,
            )
            if warning:
                warnings.append(warning)

    elif mode == "solve_from_mdot_Tt2_Pt2_M2":
        if inputs.mach_2_target <= 0 or inputs.mach_2_target >= 1:
            warnings.append("Station-2 mode B selected, but M2 target is not in the subsonic range (0,1).")
        elif mdot_2 is None or tt_2 is None or pt_2 is None:
            warnings.append("Station-2 mode B could not run because mdot_2, Tt2, or Pt2 was invalid.")
        else:
            mach_2 = inputs.mach_2_target
            try:
                a2_flow_used = area_required_for_mdot(
                    mdot_kg_s=mdot_2,
                    mach=mach_2,
                    pt_kpa=pt_2,
                    tt_k=tt_2,
                    gamma=inputs.gamma,
                    r=inputs.r_j_per_kgk,
                )
                if a2_flow_used <= 0:
                    warnings.append("Solved station-2 flow area is non-positive.")
                    a2_flow_used = None
            except Exception as e:
                warnings.append(f"Failed to solve station-2 flow area from mdot, Tt2, Pt2, and M2: {e}")
    else:
        warnings.append(
            "Invalid station2_mode. Use either "
            "'solve_from_mdot_Tt2_Pt2_A2' or 'solve_from_mdot_Tt2_Pt2_M2'."
        )

    # -------------------------------------------------------------------------
    # Derived station 2 static values
    # -------------------------------------------------------------------------
    if mach_2 is not None and tt_2 is not None and pt_2 is not None:
        try:
            t_2 = static_temperature_from_total(tt_2, mach_2, inputs.gamma)
            p_2 = static_pressure_from_total(pt_2, mach_2, inputs.gamma)

            if t_2 <= 0:
                warnings.append("Derived station-2 static temperature is non-positive.")
                t_2 = None
            if p_2 is not None and p_2 <= 0:
                warnings.append("Derived station-2 static pressure is non-positive.")
                p_2 = None

            if t_2 is not None and p_2 is not None:
                rho_2 = compute_density(p_2, t_2, inputs.r_j_per_kgk)
                a_2_sound = compute_speed_of_sound(inputs.gamma, inputs.r_j_per_kgk, t_2)
                v_2 = mach_2 * a_2_sound
        except Exception as e:
            warnings.append(f"Failed to derive station-2 static properties: {e}")

    # -------------------------------------------------------------------------
    # Ratios
    # -------------------------------------------------------------------------
    capture_to_a2flow_ideal = None
    capture_to_a2flow_working = None

    if a_capture_ideal is not None and a2_flow_used is not None and a2_flow_used > 0:
        capture_to_a2flow_ideal = a_capture_ideal / a2_flow_used

    if a_capture_working is not None and a2_flow_used is not None and a2_flow_used > 0:
        capture_to_a2flow_working = a_capture_working / a2_flow_used

    # -------------------------------------------------------------------------
    # Reference residuals
    # -------------------------------------------------------------------------
    delta_mach_2 = None if (inputs.ref_mach_2 is None or mach_2 is None) else mach_2 - inputs.ref_mach_2
    delta_t_2 = None if (inputs.ref_t_2_k is None or t_2 is None) else t_2 - inputs.ref_t_2_k
    delta_p_2 = None if (inputs.ref_p_2_kpa is None or p_2 is None) else p_2 - inputs.ref_p_2_kpa
    delta_v_2 = None if (inputs.ref_v_2_m_s is None or v_2 is None) else v_2 - inputs.ref_v_2_m_s
    delta_rho_2 = None if (inputs.ref_rho_2_kg_m3 is None or rho_2 is None) else rho_2 - inputs.ref_rho_2_kg_m3

    return IntakeSizingResults(
        rho_0_kg_m3=rho_0,
        a_0_m_s=a_0,
        v_0_m_s=v_0,
        mass_flux_0_kg_m2_s=mass_flux_0,
        tt_0_k=tt_0,
        pt_0_kpa=pt_0,
        a_capture_ideal_m2=a_capture_ideal,
        a_capture_working_m2=a_capture_working,
        d_capture_axisym_ideal_m=d_capture_axisym_ideal,
        d_capture_axisym_working_m=d_capture_axisym_working,
        area_2_input_m2=inputs.area_2_input_m2,
        area_2_input_type=inputs.area_2_input_type,
        rh_over_rt=inputs.rh_over_rt,
        rt_2_m=rt_2,
        rh_2_m=rh_2,
        d2_outer_m=d2_outer,
        a2_gross_m2=a2_gross,
        a2_annulus_m2=a2_annulus,
        a2_flow_used_m2=a2_flow_used,
        station2_mode=mode,
        mdot_2_kg_s=mdot_2,
        tt_2_k=tt_2,
        pt_2_kpa=pt_2,
        mach_2=mach_2,
        t_2_k=t_2,
        p_2_kpa=p_2,
        rho_2_kg_m3=rho_2,
        a_2_m_s=a_2_sound,
        v_2_m_s=v_2,
        capture_to_a2flow_ideal=capture_to_a2flow_ideal,
        capture_to_a2flow_working=capture_to_a2flow_working,
        two_d_examples=two_d_examples,
        delta_mach_2=delta_mach_2,
        delta_t_2_k=delta_t_2,
        delta_p_2_kpa=delta_p_2,
        delta_v_2_m_s=delta_v_2,
        delta_rho_2_kg_m3=delta_rho_2,
        warnings=warnings,
    )


def print_report(inputs: IntakeSizingInputs, results: IntakeSizingResults) -> None:
    print("=" * 104)
    print("PRELIMINARY SUPERSONIC INTAKE SIZING REPORT - V7")
    print("=" * 104)

    print("\nPURPOSE")
    print("-" * 104)
    print("First-pass intake sizing with rigid station-0 and station-2 treatment.")
    print("This version explicitly interprets the GasTurb/raw station-2 area as gross or annulus.")

    print("\nINPUTS")
    print("-" * 104)
    print(f"Station-0 Mach, M0                                           : {inputs.mach_0:.4f}")
    print(f"Cruise altitude                                              : {inputs.altitude_ft:.1f} ft")
    print(f"Station-0 static temperature, T0                             : {inputs.t_0_k:.4f} K")
    print(f"Station-0 static pressure, p0                                : {inputs.p_0_kpa:.4f} kPa")
    print(f"Per-engine actual mass flow for capture sizing               : {inputs.mdot_per_engine_kg_s:.4f} kg/s")

    print(f"\nStation-2 solve mode                                         : {inputs.station2_mode}")
    print(f"Station-2 mdot input                                         : {inputs.mdot_2_kg_s:.4f} kg/s")
    print(f"Station-2 total temperature, Tt2 input                       : {inputs.tt_2_k:.4f} K")
    print(f"Station-2 total pressure, Pt2 input                          : {inputs.pt_2_kpa:.4f} kPa")
    print(f"Station-2 raw area input                                     : {inputs.area_2_input_m2:.6f} m^2")
    print(f"Station-2 raw area interpretation                            : {inputs.area_2_input_type}")
    print(f"Station-2 hub/tip ratio, Rh/Rt                               : {inputs.rh_over_rt:.6f}")
    if inputs.station2_mode == "solve_from_mdot_Tt2_Pt2_M2":
        print(f"Station-2 Mach, M2 input                                     : {inputs.mach_2_target:.4f}")

    print(f"\nCapture area margin                                          : {100.0 * inputs.capture_area_margin_fraction:.2f} %")

    print("\nASSUMPTIONS")
    print("-" * 104)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 104)
        for i, warning in enumerate(results.warnings, start=1):
            print(f"{i:>2d}. {warning}")

    print("\nSTATION 0 STATIC")
    print("-" * 104)
    print(f"Density, rho0                                                : {fmt(results.rho_0_kg_m3, 6)} kg/m^3")
    print(f"Speed of sound, a0                                           : {fmt(results.a_0_m_s, 4)} m/s")
    print(f"Velocity, V0                                                 : {fmt(results.v_0_m_s, 4)} m/s")
    print(f"Mass flux, rho0 * V0                                         : {fmt(results.mass_flux_0_kg_m2_s, 4)} kg/(m^2 s)")

    print("\nSTATION 0 TOTAL")
    print("-" * 104)
    print(f"Total temperature, Tt0                                       : {fmt(results.tt_0_k, 4)} K")
    print(f"Total pressure, Pt0                                          : {fmt(results.pt_0_kpa, 4)} kPa")

    print("\nCAPTURE AREA PER ENGINE")
    print("-" * 104)
    print(f"Ideal capture area                                           : {fmt(results.a_capture_ideal_m2, 4)} m^2")
    print(f"Working capture area                                         : {fmt(results.a_capture_working_m2, 4)} m^2")

    print("\nAXISYMMETRIC EQUIVALENT CAPTURE SIZE PER ENGINE")
    print("-" * 104)
    print(f"Ideal equivalent capture diameter                            : {fmt(results.d_capture_axisym_ideal_m, 4)} m")
    print(f"Working equivalent capture diameter                          : {fmt(results.d_capture_axisym_working_m, 4)} m")

    print("\nSTATION 2 GEOMETRY / AREA INTERPRETATION")
    print("-" * 104)
    print(f"Raw area input                                               : {fmt(results.area_2_input_m2, 6)} m^2")
    print(f"Raw area interpreted as                                      : {results.area_2_input_type}")
    print(f"Hub/tip ratio, Rh/Rt                                         : {fmt(results.rh_over_rt, 6)}")
    print(f"Solved tip radius, Rt                                        : {fmt(results.rt_2_m, 6)} m")
    print(f"Solved hub radius, Rh                                        : {fmt(results.rh_2_m, 6)} m")
    print(f"Solved engine face outer diameter, D2                        : {fmt(results.d2_outer_m, 6)} m")
    print(f"Gross face area from solved geometry                         : {fmt(results.a2_gross_m2, 6)} m^2")
    print(f"Annulus area from solved geometry                            : {fmt(results.a2_annulus_m2, 6)} m^2")
    print(f"Flow area used in station-2 thermodynamic solve              : {fmt(results.a2_flow_used_m2, 6)} m^2")

    print("\nSTATION 2 SOLVED STATE")
    print("-" * 104)
    print(f"Station-2 mode used                                          : {results.station2_mode}")
    print(f"Mass flow, mdot2                                             : {fmt(results.mdot_2_kg_s, 4)} kg/s")
    print(f"Total temperature, Tt2                                       : {fmt(results.tt_2_k, 4)} K")
    print(f"Total pressure, Pt2                                          : {fmt(results.pt_2_kpa, 4)} kPa")
    print(f"Solved Mach, M2                                              : {fmt(results.mach_2, 6)}")
    print(f"Static temperature, T2                                       : {fmt(results.t_2_k, 4)} K")
    print(f"Static pressure, p2                                          : {fmt(results.p_2_kpa, 4)} kPa")
    print(f"Density, rho2                                                : {fmt(results.rho_2_kg_m3, 6)} kg/m^3")
    print(f"Speed of sound, a2                                           : {fmt(results.a_2_m_s, 4)} m/s")
    print(f"Velocity, V2                                                 : {fmt(results.v_2_m_s, 4)} m/s")

    print("\nAREA RATIOS")
    print("-" * 104)
    print(f"Capture / station-2 flow area (ideal)                        : {fmt(results.capture_to_a2flow_ideal, 4)}")
    print(f"Capture / station-2 flow area (working)                      : {fmt(results.capture_to_a2flow_working, 4)}")

    print("\nEXAMPLE 2D OPENINGS PER ENGINE")
    print("-" * 104)
    print(f"{'Width (m)':>12} {'Ideal height (m)':>20} {'Working height (m)':>22}")
    for width, ideal_h, working_h in results.two_d_examples:
        print(f"{width:12.4f} {fmt(ideal_h, 4):>20} {fmt(working_h, 4):>22}")

    print("\nREFERENCE CHECKS (OPTIONAL)")
    print("-" * 104)
    print(f"Delta M2                                                     : {fmt(results.delta_mach_2, 6)}")
    print(f"Delta T2 (K)                                                 : {fmt(results.delta_t_2_k, 6)}")
    print(f"Delta p2 (kPa)                                               : {fmt(results.delta_p_2_kpa, 6)}")
    print(f"Delta V2 (m/s)                                               : {fmt(results.delta_v_2_m_s, 6)}")
    print(f"Delta rho2 (kg/m^3)                                          : {fmt(results.delta_rho_2_kg_m3, 6)}")

    print("\nWHOLE AIRCRAFT TOTAL CAPTURE AREA")
    print("-" * 104)
    total_ideal = results.a_capture_ideal_m2 * inputs.engines_total if results.a_capture_ideal_m2 is not None else None
    total_working = results.a_capture_working_m2 * inputs.engines_total if results.a_capture_working_m2 is not None else None
    print(f"Total ideal capture area                                     : {fmt(total_ideal, 4)} m^2")
    print(f"Total working capture area                                   : {fmt(total_working, 4)} m^2")

    print("\nLIMITATION OF THIS VERSION")
    print("-" * 104)
    print("This script is now rigid about station definitions, area meaning, and downstream solve mode.")
    print("It still does NOT predict shock structure, throat sizing, diffuser losses,")
    print("distortion, bleed/bypass requirements, or off-design behaviour.")
    print("=" * 104)


def main() -> None:
    inputs = IntakeSizingInputs(
        # Station 0
        mach_0=1.8,
        altitude_ft=60000.0,
        t_0_k=216.65,
        p_0_kpa=7.172,

        gamma=1.4,
        r_j_per_kgk=287.05,

        engines_total=4,
        mdot_per_engine_kg_s=159.19,
        capture_area_margin_fraction=0.08,
        example_2d_widths_m=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5],

        # Station 2 mode
        station2_mode="solve_from_mdot_Tt2_Pt2_A2",

        # Station 2 independent inputs
        mdot_2_kg_s=159.19,
        tt_2_k=357.049,
        pt_2_kpa=39.169,

        # Raw area + interpretation
        area_2_input_m2=2.020917,
        area_2_input_type="gross",   # change to "annulus" if appropriate
        rh_over_rt=0.45,

        # Used only in mode B
        mach_2_target=0.748,

        # Optional reference checks
        ref_mach_2=0.748,
        ref_t_2_k=321.264,
        ref_p_2_kpa=27.0315,
        ref_v_2_m_s=268.661,
        ref_rho_2_kg_m3=0.293123,
    )

    results = compute_intake_sizing(inputs)
    print_report(inputs, results)


if __name__ == "__main__":
    main()