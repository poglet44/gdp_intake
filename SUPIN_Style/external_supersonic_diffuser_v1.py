#!/usr/bin/env python3
"""
external_supersonic_diffuser_v1.py

Purpose
-------
Script 2 in the intake design chain.

This script performs a first-pass design of the 2D external supersonic diffuser only.
It follows a station-based, literature-aligned preliminary design approach.

It imports the fixed boundary-condition / endpoint data from the station2_consistency
script and then determines:

- capture height from chosen capture width
- a 2-ramp external compression system
- EX state (end of external supersonic diffuser, just upstream of terminal shock)
- external diffuser ramp coordinates and axial length
- external shock-system total pressure recovery

What this script does NOT do
----------------------------
- apply the terminal normal shock
- compute true post-shock station 1 / NS state
- internal compression / throat section design
- subsonic diffuser design
- bifurcation duct design
- bleed / bypass / auxiliary inlet design
- viscous corrections beyond ideal oblique-shock relations

Assumptions
-----------
1. The imported Script-1 outputs are fixed design-basis inputs.
2. The inlet is treated as a 2D external supersonic diffuser in side profile.
3. The external diffuser consists of two directly connected flat ramps.
4. The script solves the state at EX, just upstream of the terminal normal shock.
5. The cowl-lip / EX end plane has a finite opening height derived from continuity.
6. Capture area is fixed and capture height is solved from chosen capture width.
7. Oblique-shock relations are ideal-gas and inviscid.
8. Equal-strength-shock criterion is approximated by matching static pressure ratios
   across the two external oblique shocks.
9. The upper cowl is horizontal in this v1 representation.
10. No explicit spillage factor is included in this v1 script.
11. Gamma is fixed at 1.4 to remain consistent with Script 1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import asin, cos, pi, sin, sqrt, tan
from typing import List, Optional, Tuple

from station2_consistency import default_inputs, compute


# --------------------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------------------

@dataclass
class ExternalDiffuserInputs:
    """
    New user-controlled inputs for the external supersonic diffuser only.
    """
    capture_width_m: float = 2.4
    MEX_target: float = 1.30
    use_working_capture_area: bool = True

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Imported Script-1 outputs are treated as fixed design-basis inputs.",
            "The external supersonic diffuser is modelled as a 2D two-ramp system.",
            "The two ramps are directly connected with no smoothing section in this v1 model.",
            "The script solves the EX state, just upstream of the terminal normal shock.",
            "Capture width is user-chosen; capture height is derived from fixed capture area.",
            "Equal-strength-shock criterion is approximated by matching the static pressure ratio across the two oblique shocks.",
            "Ideal oblique-shock relations are used; no viscous, bleed, or CFD corrections are included.",
            "This script designs the external diffuser only; throat, bifurcation, terminal shock, and subsonic diffuser are not included.",
            "Gamma is fixed at 1.4 to match Script 1.",
            "Geometry closure assumes a horizontal upper cowl from capture plane to the cowl-lip / EX end plane.",
        ]
    )


@dataclass
class ExternalDiffuserResults:
    # Imported fixed basis
    A_capture_used_m2: float
    M0: float
    T0_K: float
    p0_kPa: float
    gamma: float

    # Capture geometry
    capture_width_m: float
    capture_height_m: float

    # Solved external ramp/shock system
    delta1_deg: float
    delta2_deg: float
    beta1_deg: float
    beta2_deg: float
    M_after_ramp1: float
    MEX_actual: float

    p2_over_p1: float
    p3_over_p2: float
    external_pt_recovery: float

    # EX state (end of external supersonic diffuser, upstream of terminal shock)
    TEX_K: float
    pEX_kPa: float
    ptEX_kPa: float
    rhoEX_kg_m3: float
    aEX_m_s: float
    VEX_m_s: float
    AEX_required_m2: float
    hEX_m: float

    # Cowl-lip / end-plane geometry
    x_capture_m: float
    x_ramp_break_m: float
    x_cowl_lip_m: float

    y_capture_lower_m: float
    y_capture_upper_m: float
    y_ramp_break_m: float
    y_cowl_lip_lower_m: float
    y_cowl_lip_upper_m: float

    external_diffuser_length_m: float

    warnings: List[str]


# --------------------------------------------------------------------------------------
# Numerical helpers
# --------------------------------------------------------------------------------------

def bisect_root(func, a: float, b: float, tol: float = 1e-8, max_iter: int = 200) -> float:
    """
    Simple bisection root finder.
    Requires func(a) and func(b) to have opposite signs.
    """
    fa = func(a)
    fb = func(b)

    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        raise ValueError("Bisection interval does not bracket a root.")

    left = a
    right = b
    f_left = fa

    for _ in range(max_iter):
        mid = 0.5 * (left + right)
        f_mid = func(mid)

        if abs(f_mid) < tol or abs(right - left) < tol:
            return mid

        if f_left * f_mid <= 0.0:
            right = mid
        else:
            left = mid
            f_left = f_mid

    return 0.5 * (left + right)


def find_bracket(func, x_min: float, x_max: float, n_samples: int = 400) -> Optional[Tuple[float, float]]:
    """
    Search for a sign-changing bracket on [x_min, x_max].
    """
    prev_x = x_min
    prev_f = func(prev_x)

    for i in range(1, n_samples + 1):
        x = x_min + (x_max - x_min) * i / n_samples
        f = func(x)

        if prev_f == 0.0:
            return (prev_x, prev_x)
        if prev_f * f < 0.0:
            return (prev_x, x)

        prev_x = x
        prev_f = f

    return None


# --------------------------------------------------------------------------------------
# Compressible flow / oblique shock utilities
# --------------------------------------------------------------------------------------

def total_temperature_from_static(T_K: float, M: float, gamma: float) -> float:
    return T_K * (1.0 + 0.5 * (gamma - 1.0) * M * M)


def total_pressure_from_static(p_kPa: float, M: float, gamma: float) -> float:
    return p_kPa * (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


def static_temperature_from_total(Tt_K: float, M: float, gamma: float) -> float:
    return Tt_K / (1.0 + 0.5 * (gamma - 1.0) * M * M)


def static_pressure_from_total(pt_kPa: float, M: float, gamma: float) -> float:
    return pt_kPa / (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


def theta_beta_m_residual(beta_rad: float, M: float, delta_rad: float, gamma: float) -> float:
    """
    Residual of the theta-beta-M oblique-shock relation.

    tan(delta) = 2 cot(beta) * (M^2 sin^2(beta) - 1) / (M^2 (gamma + cos(2 beta)) + 2)
    """
    sin_b = sin(beta_rad)
    cos_2b = cos(2.0 * beta_rad)

    numerator = 2.0 * (1.0 / tan(beta_rad)) * (M * M * sin_b * sin_b - 1.0)
    denominator = M * M * (gamma + cos_2b) + 2.0
    rhs = numerator / denominator

    return tan(delta_rad) - rhs


def solve_oblique_shock_beta_weak(M: float, delta_deg: float, gamma: float) -> float:
    """
    Solve for weak-shock beta (degrees) for a given M and deflection delta.
    """
    if M <= 1.0:
        raise ValueError("Oblique shock requires upstream Mach > 1.")

    delta_rad = delta_deg * pi / 180.0
    mu = asin(1.0 / M)

    beta_min = mu + 1e-6
    beta_max = 0.5 * pi - 1e-6

    def f(beta):
        return theta_beta_m_residual(beta, M, delta_rad, gamma)

    bracket = find_bracket(f, beta_min, beta_max, n_samples=1200)
    if bracket is None:
        raise ValueError(
            f"Could not bracket a weak oblique-shock solution for M={M:.6f}, delta={delta_deg:.6f} deg."
        )

    beta_root = bisect_root(f, bracket[0], bracket[1], tol=1e-10, max_iter=300)
    return beta_root * 180.0 / pi


def normal_shock_pt_ratio_from_Mn(Mn1: float, gamma: float) -> float:
    """
    Total pressure ratio across a normal shock based on upstream normal Mach number.
    """
    if Mn1 <= 1.0:
        raise ValueError("Normal shock requires Mn1 > 1.")

    p2_over_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn1 * Mn1 - 1.0)

    Mn2_sq = (1.0 + 0.5 * (gamma - 1.0) * Mn1 * Mn1) / (gamma * Mn1 * Mn1 - 0.5 * (gamma - 1.0))
    Mn2 = sqrt(Mn2_sq)

    term_up = 1.0 + 0.5 * (gamma - 1.0) * Mn1 * Mn1
    term_down = 1.0 + 0.5 * (gamma - 1.0) * Mn2 * Mn2

    return p2_over_p1 * (term_down / term_up) ** (gamma / (gamma - 1.0))


def oblique_shock_step(M_in: float, delta_deg: float, gamma: float) -> dict:
    """
    Apply one oblique shock for a given upstream Mach and turn angle.
    """
    beta_deg = solve_oblique_shock_beta_weak(M_in, delta_deg, gamma)

    beta_rad = beta_deg * pi / 180.0
    delta_rad = delta_deg * pi / 180.0

    Mn1 = M_in * sin(beta_rad)
    if Mn1 <= 1.0:
        raise ValueError("Computed Mn1 <= 1 for oblique shock; invalid compression state.")

    p_out_over_p_in = 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn1 * Mn1 - 1.0)

    Mn2_sq = (1.0 + 0.5 * (gamma - 1.0) * Mn1 * Mn1) / (gamma * Mn1 * Mn1 - 0.5 * (gamma - 1.0))
    Mn2 = sqrt(Mn2_sq)

    M_out = Mn2 / sin(beta_rad - delta_rad)
    pt_out_over_pt_in = normal_shock_pt_ratio_from_Mn(Mn1, gamma)

    return {
        "beta_deg": beta_deg,
        "M_out": M_out,
        "p_out_over_p_in": p_out_over_p_in,
        "pt_out_over_pt_in": pt_out_over_pt_in,
        "Mn1": Mn1,
    }


# --------------------------------------------------------------------------------------
# External diffuser solve
# --------------------------------------------------------------------------------------

def compute_external_supersonic_diffuser(
    foundation_inputs,
    foundation_results,
    inputs: ExternalDiffuserInputs,
) -> ExternalDiffuserResults:
    """
    Compute the 2D 2-ramp external supersonic diffuser.
    """
    warnings: List[str] = []

    # --------------------------------------------------------------------------
    # Fixed Script-1 interface
    # --------------------------------------------------------------------------
    if inputs.use_working_capture_area:
        A_capture_used = foundation_results.A_capture_working_m2
    else:
        A_capture_used = foundation_results.A_capture_ideal_m2

    M0 = foundation_inputs.M0
    T0_K = foundation_results.T0_K
    p0_kPa = foundation_results.p0_kPa
    gamma = 1.4

    if A_capture_used is None or A_capture_used <= 0.0:
        raise ValueError("Imported capture area is invalid.")
    if M0 <= 1.0:
        raise ValueError("Freestream Mach number must be > 1 for a supersonic external diffuser.")
    if inputs.capture_width_m <= 0.0:
        raise ValueError("capture_width_m must be > 0.")
    if not (1.0 < inputs.MEX_target < M0):
        raise ValueError("MEX_target must satisfy 1 < MEX_target < M0.")

    # --------------------------------------------------------------------------
    # Capture opening
    # --------------------------------------------------------------------------
    capture_width = inputs.capture_width_m
    capture_height = A_capture_used / capture_width

    # --------------------------------------------------------------------------
    # Freestream total conditions
    # --------------------------------------------------------------------------
    Tt0_K = total_temperature_from_static(T0_K, M0, gamma)
    pt0_kPa = total_pressure_from_static(p0_kPa, M0, gamma)

    # --------------------------------------------------------------------------
    # Solve 2-ramp external compression system
    # --------------------------------------------------------------------------
    def solve_delta2_for_target_MEX(M_after_1: float, target_MEX: float) -> float:
        if target_MEX >= M_after_1:
            raise ValueError("target_MEX must be less than M_after_1.")

        def g(delta2_deg: float) -> float:
            step2 = oblique_shock_step(M_after_1, delta2_deg, gamma)
            return step2["M_out"] - target_MEX

        bracket = find_bracket(g, 0.01, 25.0, n_samples=500)
        if bracket is None:
            raise ValueError(
                f"Could not bracket delta2 for M_after_1={M_after_1:.6f}, target_MEX={target_MEX:.6f}."
            )

        return bisect_root(g, bracket[0], bracket[1], tol=1e-10, max_iter=300)

    def equal_strength_error(delta1_deg: float) -> float:
        step1 = oblique_shock_step(M0, delta1_deg, gamma)
        M_after_1 = step1["M_out"]

        if M_after_1 <= inputs.MEX_target:
            return -1e3

        delta2_deg = solve_delta2_for_target_MEX(M_after_1, inputs.MEX_target)
        step2 = oblique_shock_step(M_after_1, delta2_deg, gamma)

        return step1["p_out_over_p_in"] - step2["p_out_over_p_in"]

    bracket = find_bracket(equal_strength_error, 0.05, 20.0, n_samples=500)
    if bracket is None:
        raise ValueError(
            "Could not bracket a delta1 solution satisfying the equal-strength-shock criterion."
        )

    delta1_deg = bisect_root(equal_strength_error, bracket[0], bracket[1], tol=1e-10, max_iter=300)

    step1 = oblique_shock_step(M0, delta1_deg, gamma)
    delta2_deg = solve_delta2_for_target_MEX(step1["M_out"], inputs.MEX_target)
    step2 = oblique_shock_step(step1["M_out"], delta2_deg, gamma)

    beta1_deg = step1["beta_deg"]
    beta2_deg = step2["beta_deg"]
    M_after_ramp1 = step1["M_out"]
    MEX_actual = step2["M_out"]

    p2_over_p1 = step1["p_out_over_p_in"]
    p3_over_p2 = step2["p_out_over_p_in"]

    external_pt_recovery = step1["pt_out_over_pt_in"] * step2["pt_out_over_pt_in"]

    # --------------------------------------------------------------------------
    # EX state
    # --------------------------------------------------------------------------
    TtEX_K = Tt0_K
    ptEX_kPa = pt0_kPa * external_pt_recovery
    TEX_K = static_temperature_from_total(TtEX_K, MEX_actual, gamma)
    pEX_kPa = static_pressure_from_total(ptEX_kPa, MEX_actual, gamma)

    rhoEX_kg_m3 = (pEX_kPa * 1000.0) / (287.05 * TEX_K)
    aEX_m_s = sqrt(gamma * 287.05 * TEX_K)
    VEX_m_s = MEX_actual * aEX_m_s

    mdotEX_kg_s = foundation_inputs.mdot_capture_kg_s

    AEX_required_m2 = mdotEX_kg_s / (rhoEX_kg_m3 * VEX_m_s)
    hEX_m = AEX_required_m2 / capture_width

    # --------------------------------------------------------------------------
    # Shock-on-lip style planar geometry construction
    #
    # Coordinate convention:
    # - x = 0 at the lower capture point
    # - x increases downstream toward the cowl lip / EX end plane
    # - upper cowl is horizontal in this v1 model
    #
    # Geometry:
    # - lower capture point       = (0, 0)
    # - upper capture point       = (0, capture_height)
    # - lower cowl-lip point      = (L, capture_height - hEX)
    # - upper cowl-lip point      = (L, capture_height)
    #
    # Ramp 1:
    #   from lower capture point to ramp break
    #   slope = tan(delta1)
    #
    # Ramp 2:
    #   from ramp break to lower cowl-lip point
    #   slope = tan(delta1 + delta2)
    #
    # Shock 2:
    #   from ramp break to upper cowl-lip point
    #   global slope = tan(delta1 + delta2 + beta2)
    # --------------------------------------------------------------------------
    y_capture_lower = 0.0
    y_capture_upper = capture_height

    y_cowl_lip_upper = capture_height
    y_cowl_lip_lower = capture_height - hEX_m

    if y_cowl_lip_lower <= 0.0:
        raise ValueError(
            "Computed cowl-lip lower wall is at or below the lower capture point. "
            "This indicates the horizontal-cowl v1 geometry is not feasible for the chosen inputs."
        )

    delta1_rad = delta1_deg * pi / 180.0
    delta2_rad = delta2_deg * pi / 180.0
    beta2_rad = beta2_deg * pi / 180.0

    m1 = tan(delta1_rad)
    m2 = tan(delta1_rad + delta2_rad)
    m_shock2 = tan(delta1_rad + delta2_rad + beta2_rad)

    if m1 <= 0.0 or m2 <= 0.0 or m_shock2 <= 0.0:
        raise ValueError("Computed ramp/shock slopes are not positive.")

    if m_shock2 <= m2:
        raise ValueError(
            "Shock-2 slope is not steeper than ramp-2 slope. "
            "Shock-on-lip geometry closure is invalid."
        )

    dx_break_to_cowl_lip = hEX_m / (m_shock2 - m2)

    if dx_break_to_cowl_lip <= 0.0:
        raise ValueError("Computed break-to-cowl-lip distance is not positive.")

    y_ramp_break = y_cowl_lip_lower - m2 * dx_break_to_cowl_lip

    if y_ramp_break <= 0.0:
        raise ValueError(
            "Computed ramp-break height is not above the lower capture point. "
            "Check the chosen MEX_target and capture width."
        )

    x_ramp_break = y_ramp_break / m1
    x_cowl_lip = x_ramp_break + dx_break_to_cowl_lip
    x_capture = 0.0

    if x_ramp_break <= x_capture or x_cowl_lip <= x_ramp_break:
        raise ValueError(
            "Computed geometry is invalid: expected x_capture < x_ramp_break < x_cowl_lip."
        )

    external_diffuser_length = x_cowl_lip - x_capture

    slope_break_to_lip = (y_cowl_lip_lower - y_ramp_break) / (x_cowl_lip - x_ramp_break)
    if abs(slope_break_to_lip - m2) > 1e-8:
        warnings.append(
            "Ramp-2 geometric closure is not matched tightly; check numerical tolerance."
        )

    if abs(p2_over_p1 - p3_over_p2) > 1e-5:
        warnings.append("Equal-strength-shock criterion not matched tightly; check numerical tolerance.")

    if MEX_actual <= 1.0:
        warnings.append("Computed EX Mach is not supersonic; this is outside the intended Script-2 scope.")

    if not (0.0 < external_pt_recovery < 1.0):
        warnings.append("External shock-system total pressure recovery is outside expected bounds (0,1).")

    warnings.append(
        "Geometry closure assumes a horizontal upper cowl from capture plane to the cowl-lip / EX end plane."
    )

    return ExternalDiffuserResults(
        A_capture_used_m2=A_capture_used,
        M0=M0,
        T0_K=T0_K,
        p0_kPa=p0_kPa,
        gamma=gamma,
        capture_width_m=capture_width,
        capture_height_m=capture_height,
        delta1_deg=delta1_deg,
        delta2_deg=delta2_deg,
        beta1_deg=beta1_deg,
        beta2_deg=beta2_deg,
        M_after_ramp1=M_after_ramp1,
        MEX_actual=MEX_actual,
        p2_over_p1=p2_over_p1,
        p3_over_p2=p3_over_p2,
        external_pt_recovery=external_pt_recovery,
        TEX_K=TEX_K,
        pEX_kPa=pEX_kPa,
        ptEX_kPa=ptEX_kPa,
        rhoEX_kg_m3=rhoEX_kg_m3,
        aEX_m_s=aEX_m_s,
        VEX_m_s=VEX_m_s,
        AEX_required_m2=AEX_required_m2,
        hEX_m=hEX_m,
        x_capture_m=x_capture,
        x_ramp_break_m=x_ramp_break,
        x_cowl_lip_m=x_cowl_lip,
        y_capture_lower_m=y_capture_lower,
        y_capture_upper_m=y_capture_upper,
        y_ramp_break_m=y_ramp_break,
        y_cowl_lip_lower_m=y_cowl_lip_lower,
        y_cowl_lip_upper_m=y_cowl_lip_upper,
        external_diffuser_length_m=external_diffuser_length,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def print_external_diffuser_report(results: ExternalDiffuserResults, inputs: ExternalDiffuserInputs) -> None:
    print("=" * 100)
    print("SCRIPT 2 — EXTERNAL SUPERSONIC DIFFUSER V1")
    print("=" * 100)

    print("\nPURPOSE")
    print("-" * 100)
    print("Design the 2D external supersonic diffuser only, using fixed Script-1 boundary conditions.")
    print("This script solves the EX state and the cowl-lip / EX end-plane geometry.")
    print("It does not yet apply the terminal normal shock to reach the true post-shock station-1 / NS state.")

    print("\nASSUMPTIONS")
    print("-" * 100)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED FIXED DESIGN BASIS")
    print("-" * 100)
    print(f"Freestream Mach, M0                             : {results.M0:.6f}")
    print(f"Freestream static temperature, T0               : {results.T0_K:.6f} K")
    print(f"Freestream static pressure, p0                  : {results.p0_kPa:.6f} kPa")
    print(f"Gamma                                           : {results.gamma:.6f}")
    print(f"Capture area used                               : {results.A_capture_used_m2:.6f} m^2")

    print("\nUSER-SET INPUTS")
    print("-" * 100)
    print(f"Capture width                                   : {inputs.capture_width_m:.6f} m")
    print(f"Target EX Mach, MEX_target                      : {inputs.MEX_target:.6f}")
    print(f"Use working capture area                        : {inputs.use_working_capture_area}")

    print("\nCAPTURE GEOMETRY")
    print("-" * 100)
    print(f"Capture width                                   : {results.capture_width_m:.6f} m")
    print(f"Capture height                                  : {results.capture_height_m:.6f} m")

    print("\nSOLVED 2-RAMP EXTERNAL COMPRESSION SYSTEM")
    print("-" * 100)
    print(f"Ramp 1 deflection, delta1                       : {results.delta1_deg:.6f} deg")
    print(f"Ramp 2 deflection, delta2                       : {results.delta2_deg:.6f} deg")
    print(f"Shock 1 angle, beta1                            : {results.beta1_deg:.6f} deg")
    print(f"Shock 2 angle, beta2                            : {results.beta2_deg:.6f} deg")
    print(f"Mach after ramp 1                               : {results.M_after_ramp1:.6f}")
    print(f"EX Mach, MEX_actual                             : {results.MEX_actual:.6f}")
    print(f"Shock 1 static pressure ratio                   : {results.p2_over_p1:.6f}")
    print(f"Shock 2 static pressure ratio                   : {results.p3_over_p2:.6f}")
    print(f"External shock-system pt recovery               : {results.external_pt_recovery:.6f}")

    print("\nEX STATE (JUST UPSTREAM OF TERMINAL NORMAL SHOCK)")
    print("-" * 100)
    print(f"EX static temperature                           : {results.TEX_K:.6f} K")
    print(f"EX static pressure                              : {results.pEX_kPa:.6f} kPa")
    print(f"EX total pressure                               : {results.ptEX_kPa:.6f} kPa")
    print(f"EX density                                      : {results.rhoEX_kg_m3:.6f} kg/m^3")
    print(f"EX speed of sound                               : {results.aEX_m_s:.6f} m/s")
    print(f"EX velocity                                     : {results.VEX_m_s:.6f} m/s")
    print(f"EX required flow area                           : {results.AEX_required_m2:.6f} m^2")
    print(f"EX / cowl-lip opening height                    : {results.hEX_m:.6f} m")

    print("\nPLANAR SIDE-VIEW GEOMETRY")
    print("-" * 100)
    print(f"Capture x                                       : {results.x_capture_m:.6f} m")
    print(f"Ramp-break x                                    : {results.x_ramp_break_m:.6f} m")
    print(f"Cowl-lip / EX end-plane x                       : {results.x_cowl_lip_m:.6f} m")
    print(f"Capture lower y                                 : {results.y_capture_lower_m:.6f} m")
    print(f"Capture upper y                                 : {results.y_capture_upper_m:.6f} m")
    print(f"Ramp-break y                                    : {results.y_ramp_break_m:.6f} m")
    print(f"Cowl-lip lower y                                : {results.y_cowl_lip_lower_m:.6f} m")
    print(f"Cowl-lip upper y                                : {results.y_cowl_lip_upper_m:.6f} m")
    print(f"External diffuser axial length                  : {results.external_diffuser_length_m:.6f} m")

    print("\nLIMITATION OF THIS STEP")
    print("-" * 100)
    print("This script defines only the 2D external supersonic diffuser up to EX / the cowl-lip end plane.")
    print("It does not yet apply the terminal normal shock, and it does not define throat geometry,")
    print("bifurcation geometry, or the subsonic diffuser.")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 100)
        for i, warning in enumerate(results.warnings, start=1):
            print(f"{i:>2d}. {warning}")

    print("=" * 100)


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------

def main() -> None:
    foundation_inputs = default_inputs()
    foundation_results = compute(foundation_inputs)

    inputs = ExternalDiffuserInputs()

    results = compute_external_supersonic_diffuser(foundation_inputs, foundation_results, inputs)
    print_external_diffuser_report(results, inputs)


if __name__ == "__main__":
    main()