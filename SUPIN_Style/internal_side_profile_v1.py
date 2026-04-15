#!/usr/bin/env python3
"""
post_lip_side_profile_check_v1.py

Purpose
-------
Script 3 in the intake design chain.

This script performs a very simple post-lip internal side-profile continuation check.

It does NOT solve:
- terminal normal shock
- post-shock station state
- throat section
- TX / SD stations
- bifurcation
- subsonic diffuser

It only asks:

"Starting from the cowl-lip / EX end-plane geometry from Script 2, if the upper wall
comes down and the lower wall rises over a short downstream distance, does the internal
opening remain sensible?"

This is a geometry sanity-check script only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi, tan
from typing import List

from station2_consistency import default_inputs, compute
from external_supersonic_diffuser_v1 import (
    ExternalDiffuserInputs,
    compute_external_supersonic_diffuser,
)


# --------------------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------------------

@dataclass
class PostLipCheckInputs:
    """
    Simple geometry factors for a first post-lip continuation check.
    """
    theta_upper_internal_deg: float = 6.0
    theta_lower_internal_deg: float = 3.0
    x_check_m: float = 1.0

    assumptions: List[str] = field(
        default_factory=lambda: [
            "This is a side-view geometry check only.",
            "The cowl-lip / EX end-plane geometry from Script 2 is treated as fixed.",
            "The upper internal wall slopes downward after the lip with theta_upper_internal_deg.",
            "The lower internal wall slopes upward after the lip with theta_lower_internal_deg.",
            "Width is constant and imported from Script 2.",
            "Area is computed as width x opening height.",
            "This script does not solve terminal shock, throat, bifurcation, or diffuser physics.",
        ]
    )


@dataclass
class PostLipCheckResults:
    # Imported basis
    capture_width_m: float
    x_lip_m: float
    y_lip_lower_m: float
    y_lip_upper_m: float
    h_lip_m: float
    A_lip_m2: float

    # Check plane
    x_check_m: float
    y_check_lower_m: float
    y_check_upper_m: float
    h_check_m: float
    A_check_m2: float

    # Trends
    delta_h_m: float
    delta_A_m2: float
    percent_height_change: float
    percent_area_change: float
    interpretation: str

    # Inputs echoed
    theta_upper_internal_deg: float
    theta_lower_internal_deg: float

    warnings: List[str]


# --------------------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------------------

def compute_post_lip_check(external_results, inputs: PostLipCheckInputs) -> PostLipCheckResults:
    """
    Compute a simple post-lip side-profile continuation check.
    """
    warnings: List[str] = []

    if inputs.x_check_m <= 0.0:
        raise ValueError("x_check_m must be > 0.")
    if inputs.theta_upper_internal_deg < 0.0:
        raise ValueError("theta_upper_internal_deg must be >= 0.")
    if inputs.theta_lower_internal_deg < 0.0:
        raise ValueError("theta_lower_internal_deg must be >= 0.")

    width = external_results.capture_width_m

    # Imported lip geometry from Script 2
    x_lip = external_results.x_cowl_lip_m
    y_lip_lower = external_results.y_cowl_lip_lower_m
    y_lip_upper = external_results.y_cowl_lip_upper_m

    h_lip = y_lip_upper - y_lip_lower
    if h_lip <= 0.0:
        raise ValueError("Imported lip opening height must be positive.")

    A_lip = width * h_lip

    # Slopes
    m_upper = tan(inputs.theta_upper_internal_deg * pi / 180.0)
    m_lower = tan(inputs.theta_lower_internal_deg * pi / 180.0)

    # Check plane downstream of lip
    x_check = x_lip + inputs.x_check_m

    y_check_upper = y_lip_upper - m_upper * inputs.x_check_m
    y_check_lower = y_lip_lower + m_lower * inputs.x_check_m

    h_check = y_check_upper - y_check_lower
    A_check = width * h_check

    delta_h = h_check - h_lip
    delta_A = A_check - A_lip

    percent_height_change = 100.0 * delta_h / h_lip
    percent_area_change = 100.0 * delta_A / A_lip

    # Interpretation
    if h_check <= 0.0:
        interpretation = "invalid"
        warnings.append("Opening becomes non-positive at the check plane.")
    else:
        reduction_fraction = (h_lip - h_check) / h_lip

        if reduction_fraction < 0.0:
            interpretation = "expanding"
        elif reduction_fraction < 0.15:
            interpretation = "mild contraction"
        elif reduction_fraction < 0.30:
            interpretation = "moderate contraction"
        else:
            interpretation = "aggressive contraction"

    # Checks / warnings
    if y_check_upper <= y_check_lower:
        warnings.append("Upper and lower walls cross by the check plane.")
    if h_check <= 0.0:
        warnings.append("Check-plane opening height is non-positive.")
    if interpretation == "aggressive contraction":
        warnings.append("Post-lip continuation looks strongly contracting and needs critique.")
    if interpretation == "expanding":
        warnings.append("Post-lip continuation expands rather than contracts.")
    if inputs.theta_upper_internal_deg == 0.0 and inputs.theta_lower_internal_deg == 0.0:
        warnings.append("Both walls are flat; this is only a neutral reference case.")

    return PostLipCheckResults(
        capture_width_m=width,
        x_lip_m=x_lip,
        y_lip_lower_m=y_lip_lower,
        y_lip_upper_m=y_lip_upper,
        h_lip_m=h_lip,
        A_lip_m2=A_lip,
        x_check_m=x_check,
        y_check_lower_m=y_check_lower,
        y_check_upper_m=y_check_upper,
        h_check_m=h_check,
        A_check_m2=A_check,
        delta_h_m=delta_h,
        delta_A_m2=delta_A,
        percent_height_change=percent_height_change,
        percent_area_change=percent_area_change,
        interpretation=interpretation,
        theta_upper_internal_deg=inputs.theta_upper_internal_deg,
        theta_lower_internal_deg=inputs.theta_lower_internal_deg,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def print_post_lip_check_report(results: PostLipCheckResults, inputs: PostLipCheckInputs) -> None:
    print("=" * 100)
    print("SCRIPT 3 — POST-LIP SIDE-PROFILE CONTINUATION CHECK V1")
    print("=" * 100)

    print("\nPURPOSE")
    print("-" * 100)
    print("Check whether the internal side profile can continue sensibly just after the cowl lip.")
    print("This is a geometry sanity check only.")

    print("\nASSUMPTIONS")
    print("-" * 100)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED LIP / EX END-PLANE BASIS")
    print("-" * 100)
    print(f"Width                                           : {results.capture_width_m:.6f} m")
    print(f"Lip x                                           : {results.x_lip_m:.6f} m")
    print(f"Lip lower y                                     : {results.y_lip_lower_m:.6f} m")
    print(f"Lip upper y                                     : {results.y_lip_upper_m:.6f} m")
    print(f"Lip opening height                              : {results.h_lip_m:.6f} m")
    print(f"Lip area                                        : {results.A_lip_m2:.6f} m^2")

    print("\nUSER-SET CHECK FACTORS")
    print("-" * 100)
    print(f"Upper-wall internal angle                       : {results.theta_upper_internal_deg:.6f} deg")
    print(f"Lower-wall internal angle                       : {results.theta_lower_internal_deg:.6f} deg")
    print(f"Check distance downstream from lip              : {inputs.x_check_m:.6f} m")

    print("\nCHECK PLANE")
    print("-" * 100)
    print(f"Check-plane x                                   : {results.x_check_m:.6f} m")
    print(f"Check-plane lower y                             : {results.y_check_lower_m:.6f} m")
    print(f"Check-plane upper y                             : {results.y_check_upper_m:.6f} m")
    print(f"Check-plane opening height                      : {results.h_check_m:.6f} m")
    print(f"Check-plane area                                : {results.A_check_m2:.6f} m^2")

    print("\nCHANGE FROM LIP TO CHECK PLANE")
    print("-" * 100)
    print(f"Height change                                   : {results.delta_h_m:.6f} m")
    print(f"Area change                                     : {results.delta_A_m2:.6f} m^2")
    print(f"Percent height change                           : {results.percent_height_change:.6f} %")
    print(f"Percent area change                             : {results.percent_area_change:.6f} %")
    print(f"Interpretation                                  : {results.interpretation}")

    print("\nHOW TO TAKE THIS RESULT")
    print("-" * 100)
    print("This script only tells you whether the internal opening just after the lip")
    print("looks mild, moderate, aggressive, expanding, or invalid.")
    print("It does not yet tell you whether the full intake is aerodynamically good.")

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
    # Script 1 basis
    foundation_inputs = default_inputs()
    foundation_results = compute(foundation_inputs)

    # Script 2 basis
    external_inputs = ExternalDiffuserInputs()
    external_results = compute_external_supersonic_diffuser(
        foundation_inputs,
        foundation_results,
        external_inputs,
    )

    # Script 3 inputs
    post_lip_inputs = PostLipCheckInputs()

    post_lip_results = compute_post_lip_check(external_results, post_lip_inputs)
    print_post_lip_check_report(post_lip_results, post_lip_inputs)


if __name__ == "__main__":
    main()