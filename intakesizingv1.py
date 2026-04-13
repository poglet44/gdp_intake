#!/usr/bin/env python3
"""
prelim_intake_sizing_v2.py

Preliminary supersonic intake sizing tool
=========================================

Purpose
-------
This script performs first-pass sizing for a supersonic aircraft intake.
It currently covers:
1. Freestream properties at the design point
2. Capture area sizing using actual engine mass flow
3. Engine/fan-face annulus geometry
4. Fan-face flow state using a target fan-face Mach number
5. Area ratios between capture plane and engine/fan face
6. Equivalent axisymmetric and 2D realised opening sizes

This is still a low-order concept tool, not a full intake design solver.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional


# =============================================================================
# ASSUMPTIONS USED IN THIS VERSION
# =============================================================================
# 1. Single design-point sizing only.
#    The intake is sized at one chosen operating condition, not across the whole mission.
#
# 2. Uniform freestream at the intake capture plane.
#    No aircraft forebody, wing, local installation, or AoA effects are included yet.
#
# 3. Capture area is based on continuity using ACTUAL engine mass flow:
#       A_capture = mdot / (rho_inf * V_inf)
#
# 4. Corrected / referred mass flow (e.g. WRstd) is NOT used directly for capture sizing.
#    If corrected flow is used in future, it must first be converted back to actual mass flow.
#
# 5. Fan/engine face is represented using:
#    - outer diameter
#    - hub/tip ratio
#    - target fan-face Mach number
#
# 6. The fan-face annulus area is treated as the available flow area at the downstream boundary.
#
# 7. No explicit modelling yet of:
#    - shock system
#    - total pressure recovery
#    - throat sizing
#    - subsonic diffuser losses
#    - bleed / bypass
#    - distortion
#    - 3D bifurcation losses
#    - off-design operability
#
# 8. A simple user-defined sizing margin is applied to move from ideal to working capture area.
#
# 9. For the same required mass flow and freestream condition, different inlet families
#    are assumed to need approximately the same capture area. Therefore, at this stage,
#    2D bifurcated and axisymmetric concepts differ mainly in how that area is realised.
#
# 10. Perfect-gas behaviour is assumed with constant gamma and R.
#
# 11. Station-2 / fan-face static conditions are estimated from the target fan-face Mach number
#     and freestream total conditions. This is only a crude first-pass approximation, because
#     real inlet total pressure losses are not yet modelled.
# =============================================================================


@dataclass
class IntakeSizingInputs:
    # -------------------------
    # Design point
    # -------------------------
    mach: float = 1.8
    altitude_ft: float = 60000.0

    # Freestream static conditions at design point
    t_inf_k: float = 216.65
    p_inf_kpa: float = 7.172

    # Perfect gas properties
    gamma: float = 1.4
    r_j_per_kgk: float = 287.05

    # -------------------------
    # Engine / intake inputs
    # -------------------------
    engines_total: int = 4
    mdot_per_engine_kg_s: float = 159.19

    # Downstream boundary
    fan_face_diameter_m: float = 1.516185702
    hub_tip_ratio: float = 0.45
    fan_face_mach_target: float = 0.6

    # Simple practical sizing margin
    capture_area_margin_fraction: float = 0.08

    # Example 2D capture widths to inspect
    example_2d_widths_m: List[float] = field(
        default_factory=lambda: [2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
    )

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Single design-point preliminary sizing only.",
            "Uniform freestream assumed at the intake capture plane.",
            "Capture area computed using actual mass flow only, not corrected/referred flow.",
            "No shock, recovery, throat, diffuser, bleed, bypass, or distortion modelling in this version.",
            "Fan-face annulus defined from outer diameter and hub/tip ratio.",
            "Fan-face state estimated using target fan-face Mach and freestream total conditions as a crude first pass.",
            "A simple user-defined margin is applied to move from ideal to working capture area.",
            "2D bifurcated and axisymmetric concepts are assumed to need approximately the same capture area for the same mdot and freestream state.",
            "Perfect-gas model with constant gamma and gas constant R.",
        ]
    )


@dataclass
class IntakeSizingResults:
    # Freestream static
    rho_inf_kg_m3: float
    a_inf_m_s: float
    v_inf_m_s: float
    mass_flux_kg_m2_s: float

    # Freestream total
    tt_inf_k: float
    pt_inf_kpa: float

    # Capture sizing
    capture_area_ideal_m2: float
    capture_area_working_m2: float
    axisymmetric_capture_diameter_ideal_m: float
    axisymmetric_capture_diameter_working_m: float

    # Fan face geometry
    hub_diameter_m: float
    face_gross_area_m2: float
    face_annulus_area_m2: float

    # Fan face first-pass state
    t_face_static_k: float
    p_face_static_kpa: float
    a_face_m_s: float
    v_face_m_s: float
    rho_face_kg_m3: float

    # Flow check at face
    mdot_face_check_kg_s: float

    # Area ratios
    capture_to_gross_face_ideal: float
    capture_to_gross_face_working: float
    capture_to_annulus_ideal: float
    capture_to_annulus_working: float

    # 2D examples: (width, ideal height, working height)
    two_d_examples: List[tuple[float, float, float]]


def circle_area(diameter_m: float) -> float:
    return math.pi * diameter_m**2 / 4.0


def equivalent_diameter_from_area(area_m2: float) -> float:
    return math.sqrt(4.0 * area_m2 / math.pi)


def compute_density(p_kpa: float, t_k: float, r: float) -> float:
    return (p_kpa * 1000.0) / (r * t_k)


def compute_speed_of_sound(gamma: float, r: float, t_k: float) -> float:
    return math.sqrt(gamma * r * t_k)


def total_temperature(t_static: float, mach: float, gamma: float) -> float:
    return t_static * (1.0 + 0.5 * (gamma - 1.0) * mach**2)


def total_pressure(p_static: float, mach: float, gamma: float) -> float:
    return p_static * (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (gamma / (gamma - 1.0))


def static_temperature_from_total(tt: float, mach: float, gamma: float) -> float:
    return tt / (1.0 + 0.5 * (gamma - 1.0) * mach**2)


def static_pressure_from_total(pt: float, mach: float, gamma: float) -> float:
    return pt / (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (gamma / (gamma - 1.0))


def fmt(value: Optional[float], ndp: int = 4) -> str:
    if value is None:
        return "N/A"
    return f"{value:.{ndp}f}"


def compute_intake_sizing(inputs: IntakeSizingInputs) -> IntakeSizingResults:
    # -------------------------------------------------------------------------
    # Freestream static properties
    # -------------------------------------------------------------------------
    rho_inf = compute_density(inputs.p_inf_kpa, inputs.t_inf_k, inputs.r_j_per_kgk)
    a_inf = compute_speed_of_sound(inputs.gamma, inputs.r_j_per_kgk, inputs.t_inf_k)
    v_inf = inputs.mach * a_inf
    mass_flux = rho_inf * v_inf

    # -------------------------------------------------------------------------
    # Freestream total properties
    # -------------------------------------------------------------------------
    tt_inf = total_temperature(inputs.t_inf_k, inputs.mach, inputs.gamma)
    pt_inf = total_pressure(inputs.p_inf_kpa, inputs.mach, inputs.gamma)

    # -------------------------------------------------------------------------
    # Capture sizing per engine
    # -------------------------------------------------------------------------
    capture_area_ideal = inputs.mdot_per_engine_kg_s / mass_flux
    capture_area_working = capture_area_ideal * (1.0 + inputs.capture_area_margin_fraction)

    axisym_d_ideal = equivalent_diameter_from_area(capture_area_ideal)
    axisym_d_working = equivalent_diameter_from_area(capture_area_working)

    # -------------------------------------------------------------------------
    # Fan-face geometry
    # -------------------------------------------------------------------------
    hub_diameter = inputs.hub_tip_ratio * inputs.fan_face_diameter_m
    face_gross_area = circle_area(inputs.fan_face_diameter_m)
    hub_area = circle_area(hub_diameter)
    face_annulus_area = face_gross_area - hub_area

    # -------------------------------------------------------------------------
    # First-pass fan-face state
    #
    # IMPORTANT:
    # This uses freestream total conditions as a crude first-pass stand-in for
    # face total conditions. Real inlet total pressure loss is not yet modelled.
    # -------------------------------------------------------------------------
    t_face_static = static_temperature_from_total(
        tt_inf, inputs.fan_face_mach_target, inputs.gamma
    )
    p_face_static = static_pressure_from_total(
        pt_inf, inputs.fan_face_mach_target, inputs.gamma
    )
    a_face = compute_speed_of_sound(inputs.gamma, inputs.r_j_per_kgk, t_face_static)
    v_face = inputs.fan_face_mach_target * a_face
    rho_face = compute_density(p_face_static, t_face_static, inputs.r_j_per_kgk)

    # Check the mass flow that such a face state would pass through the annulus
    mdot_face_check = rho_face * v_face * face_annulus_area

    # -------------------------------------------------------------------------
    # Area ratios
    # -------------------------------------------------------------------------
    capture_to_gross_face_ideal = capture_area_ideal / face_gross_area
    capture_to_gross_face_working = capture_area_working / face_gross_area
    capture_to_annulus_ideal = capture_area_ideal / face_annulus_area
    capture_to_annulus_working = capture_area_working / face_annulus_area

    # -------------------------------------------------------------------------
    # Example 2D realised capture openings
    # -------------------------------------------------------------------------
    two_d_examples = []
    for width in inputs.example_2d_widths_m:
        ideal_height = capture_area_ideal / width
        working_height = capture_area_working / width
        two_d_examples.append((width, ideal_height, working_height))

    return IntakeSizingResults(
        rho_inf_kg_m3=rho_inf,
        a_inf_m_s=a_inf,
        v_inf_m_s=v_inf,
        mass_flux_kg_m2_s=mass_flux,
        tt_inf_k=tt_inf,
        pt_inf_kpa=pt_inf,
        capture_area_ideal_m2=capture_area_ideal,
        capture_area_working_m2=capture_area_working,
        axisymmetric_capture_diameter_ideal_m=axisym_d_ideal,
        axisymmetric_capture_diameter_working_m=axisym_d_working,
        hub_diameter_m=hub_diameter,
        face_gross_area_m2=face_gross_area,
        face_annulus_area_m2=face_annulus_area,
        t_face_static_k=t_face_static,
        p_face_static_kpa=p_face_static,
        a_face_m_s=a_face,
        v_face_m_s=v_face,
        rho_face_kg_m3=rho_face,
        mdot_face_check_kg_s=mdot_face_check,
        capture_to_gross_face_ideal=capture_to_gross_face_ideal,
        capture_to_gross_face_working=capture_to_gross_face_working,
        capture_to_annulus_ideal=capture_to_annulus_ideal,
        capture_to_annulus_working=capture_to_annulus_working,
        two_d_examples=two_d_examples,
    )


def print_report(inputs: IntakeSizingInputs, results: IntakeSizingResults) -> None:
    print("=" * 88)
    print("PRELIMINARY SUPERSONIC INTAKE SIZING REPORT - V2")
    print("=" * 88)

    print("\nPURPOSE")
    print("-" * 88)
    print("First-pass intake sizing from capture plane to fan/engine face for early concept work.")

    print("\nINPUTS")
    print("-" * 88)
    print(f"Cruise Mach                                      : {inputs.mach:.4f}")
    print(f"Cruise altitude                                  : {inputs.altitude_ft:.1f} ft")
    print(f"Freestream temperature                           : {inputs.t_inf_k:.4f} K")
    print(f"Freestream pressure                              : {inputs.p_inf_kpa:.4f} kPa")
    print(f"Per-engine actual mass flow                      : {inputs.mdot_per_engine_kg_s:.4f} kg/s")
    print(f"Total number of engines                          : {inputs.engines_total:d}")
    print(f"Fan/engine face outer diameter                   : {inputs.fan_face_diameter_m:.6f} m")
    print(f"Hub/tip ratio                                    : {inputs.hub_tip_ratio:.4f}")
    print(f"Fan face Mach target                             : {inputs.fan_face_mach_target:.4f}")
    print(f"Capture area margin                              : {100.0 * inputs.capture_area_margin_fraction:.2f} %")

    print("\nASSUMPTIONS")
    print("-" * 88)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nFREESTREAM STATIC")
    print("-" * 88)
    print(f"Density, rho_inf                                 : {results.rho_inf_kg_m3:.6f} kg/m^3")
    print(f"Speed of sound, a_inf                            : {results.a_inf_m_s:.4f} m/s")
    print(f"Flight speed, V_inf                              : {results.v_inf_m_s:.4f} m/s")
    print(f"Mass flux, rho_inf * V_inf                       : {results.mass_flux_kg_m2_s:.4f} kg/(m^2 s)")

    print("\nFREESTREAM TOTAL")
    print("-" * 88)
    print(f"Total temperature, Tt_inf                        : {results.tt_inf_k:.4f} K")
    print(f"Total pressure, Pt_inf                           : {results.pt_inf_kpa:.4f} kPa")

    print("\nCAPTURE AREA PER ENGINE")
    print("-" * 88)
    print(f"Ideal capture area                               : {results.capture_area_ideal_m2:.4f} m^2")
    print(f"Working capture area                             : {results.capture_area_working_m2:.4f} m^2")

    print("\nAXISYMMETRIC EQUIVALENT CAPTURE SIZE PER ENGINE")
    print("-" * 88)
    print(f"Ideal equivalent capture diameter                : {results.axisymmetric_capture_diameter_ideal_m:.4f} m")
    print(f"Working equivalent capture diameter              : {results.axisymmetric_capture_diameter_working_m:.4f} m")

    print("\nFAN / ENGINE FACE GEOMETRY")
    print("-" * 88)
    print(f"Hub diameter                                     : {results.hub_diameter_m:.4f} m")
    print(f"Gross face area                                  : {results.face_gross_area_m2:.4f} m^2")
    print(f"Annulus face area                                : {results.face_annulus_area_m2:.4f} m^2")

    print("\nFIRST-PASS FAN/FACE FLOW STATE")
    print("-" * 88)
    print("WARNING: this is based on freestream total conditions with no inlet losses yet.")
    print(f"Target face Mach                                 : {inputs.fan_face_mach_target:.4f}")
    print(f"Static temperature at face                       : {results.t_face_static_k:.4f} K")
    print(f"Static pressure at face                          : {results.p_face_static_kpa:.4f} kPa")
    print(f"Speed of sound at face                           : {results.a_face_m_s:.4f} m/s")
    print(f"Velocity at face                                 : {results.v_face_m_s:.4f} m/s")
    print(f"Density at face                                  : {results.rho_face_kg_m3:.6f} kg/m^3")
    print(f"Mass flow passed by annulus at this state        : {results.mdot_face_check_kg_s:.4f} kg/s")

    print("\nAREA RATIOS")
    print("-" * 88)
    print(f"Capture / gross face area (ideal)                : {results.capture_to_gross_face_ideal:.4f}")
    print(f"Capture / gross face area (working)              : {results.capture_to_gross_face_working:.4f}")
    print(f"Capture / annulus area (ideal)                   : {results.capture_to_annulus_ideal:.4f}")
    print(f"Capture / annulus area (working)                 : {results.capture_to_annulus_working:.4f}")

    print("\nEXAMPLE 2D OPENINGS PER ENGINE")
    print("-" * 88)
    print(f"{'Width (m)':>12} {'Ideal height (m)':>20} {'Working height (m)':>22}")
    for width, ideal_h, working_h in results.two_d_examples:
        print(f"{width:12.4f} {ideal_h:20.4f} {working_h:22.4f}")

    print("\nWHOLE AIRCRAFT TOTAL CAPTURE AREA")
    print("-" * 88)
    total_ideal = results.capture_area_ideal_m2 * inputs.engines_total
    total_working = results.capture_area_working_m2 * inputs.engines_total
    print(f"Total ideal capture area                         : {total_ideal:.4f} m^2")
    print(f"Total working capture area                       : {total_working:.4f} m^2")

    print("\nLIMITATION OF THIS VERSION")
    print("-" * 88)
    print("This script still does NOT predict shock structure, pressure recovery,")
    print("throat area, subsonic diffuser losses, distortion, or off-design behaviour.")
    print("It is now suitable for first-pass sizing from capture plane to fan/engine face only.")

    print("=" * 88)


def main() -> None:
    inputs = IntakeSizingInputs(
        mach=1.8,
        altitude_ft=60000.0,
        t_inf_k=216.65,
        p_inf_kpa=7.172,
        gamma=1.4,
        r_j_per_kgk=287.05,
        engines_total=4,
        mdot_per_engine_kg_s=159.19,
        fan_face_diameter_m=1.516185702,
        hub_tip_ratio=0.45,
        fan_face_mach_target=0.6,
        capture_area_margin_fraction=0.08,
        example_2d_widths_m=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5],
    )

    results = compute_intake_sizing(inputs)
    print_report(inputs, results)


if __name__ == "__main__":
    main()