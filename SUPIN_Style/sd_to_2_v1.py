#!/usr/bin/env python3
"""
sd_to_2_v2.py

Purpose
-------
Script 6 in the intake design chain.

This script solves the subsonic diffuser from station SD to station 2 with a
shape-transition-aware first-pass model.

Main path
---------
0 -> EX -> 1/NS -> TH -> SD -> 2

This script computes:
- station-2 equivalent side-view target geometry
- SD -> 2 area / height / width progression
- first-pass subsonic diffuser length using a shape-transition-aware rule
- first-pass end coordinates

This script does NOT compute:
- detailed blended wall surfaces
- detailed super-ellipse definitions
- bifurcation integration
- CFD mesh geometry

It is a station/area/width/height/length script only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi, sqrt, tan
from typing import List

from station2_consistency import default_inputs, compute
from external_supersonic_diffuser_v1 import (
    ExternalDiffuserInputs,
    compute_external_supersonic_diffuser,
)
from terminal_normal_shock_v1 import compute_terminal_normal_shock
from ns_th_sd_v1 import compute_ns_th_sd, NSTHSDInputs


# --------------------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------------------

@dataclass
class SDTo2Inputs:
    """
    New design factors for SD -> 2.
    """
    equivalent_conical_angle_deg: float = 3.0
    equivalent_width_transition_deg: float = 6.0
    equivalent_height_transition_deg: float = 6.0
    width_mode: str = "linear_to_D2"   # "constant" or "linear_to_D2"

    assumptions: List[str] = field(
        default_factory=lambda: [
            "This script solves SD -> 2 in area/height/width/length terms only.",
            "The subsonic diffuser is the section that carries the rectangular-to-circular/annular transition for a 2D inlet.",
            "A first-pass diffuser length is computed using a shape-transition-aware rule.",
            "The controlling SD->2 length is taken as the maximum of the area-based, width-based, and height-based length estimates.",
            "The first-pass geometry is represented using equivalent widths and heights, not full blended surfaces.",
            "This script does not yet generate the detailed blended subsonic diffuser surfaces needed for CFD.",
        ]
    )


@dataclass
class SDTo2Results:
    # Imported SD basis
    ASD_m2: float
    hSD_m: float
    wSD_m: float

    # Imported station-2 basis
    A2_m2: float
    D2_m: float
    h2_eq_m: float
    w2_eq_m: float

    # Candidate lengths
    L_area_based_m: float
    L_width_based_m: float
    L_height_based_m: float

    # Chosen length
    L_SD_to_2_m: float
    controlling_mode: str

    # Coordinates
    x_SD_m: float
    x_2_m: float

    # Ratios
    A2_over_ASD: float
    h2_over_hSD: float
    w2_over_wSD: float
    percent_height_change: float
    percent_width_change: float

    interpretation_area: str
    interpretation_width: str
    interpretation_height: str

    warnings: List[str]


# --------------------------------------------------------------------------------------
# Core helpers
# --------------------------------------------------------------------------------------

def interpret_change(end_val: float, start_val: float, label: str) -> str:
    frac = (end_val - start_val) / start_val
    if frac > 0.15:
        return f"aggressive {label} increase"
    if frac > 0.03:
        return f"moderate {label} increase"
    if frac >= -0.03:
        return f"near-constant {label}"
    if frac >= -0.15:
        return f"moderate {label} decrease"
    return f"aggressive {label} decrease"


def compute_area_based_length(
    ASD_m2: float,
    A2_m2: float,
    equivalent_conical_angle_deg: float,
) -> float:
    """
    First-pass diffuser length from equivalent conical-angle rule.

    Uses:
        r_eq = sqrt(A/pi)
        L = |r2 - r1| / tan(theta)
    """
    if equivalent_conical_angle_deg <= 0.0:
        raise ValueError("equivalent_conical_angle_deg must be > 0.")

    rSD_eq = sqrt(ASD_m2 / pi)
    r2_eq = sqrt(A2_m2 / pi)

    delta_r = abs(r2_eq - rSD_eq)
    if delta_r == 0.0:
        return 0.0

    return delta_r / tan(equivalent_conical_angle_deg * pi / 180.0)


def compute_linear_dimension_length(
    start_val: float,
    end_val: float,
    equivalent_angle_deg: float,
) -> float:
    """
    Convert a width or height change into a first-pass axial length using
    an equivalent transition angle.
    """
    if equivalent_angle_deg <= 0.0:
        raise ValueError("equivalent_angle_deg must be > 0.")

    delta_val = abs(end_val - start_val)
    if delta_val == 0.0:
        return 0.0

    return delta_val / tan(equivalent_angle_deg * pi / 180.0)


# --------------------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------------------

def compute_sd_to_2(
    ns_th_sd_results,
    foundation_results,
    inputs: SDTo2Inputs,
) -> SDTo2Results:
    """
    Solve SD -> 2 in area/height/width/length terms with a
    shape-transition-aware length model.
    """
    warnings: List[str] = []

    ASD_m2 = ns_th_sd_results.ASD_required_m2
    hSD_m = ns_th_sd_results.hSD_m
    wSD_m = ns_th_sd_results.width_m

    A2_m2 = foundation_results.A2_annulus_req_m2
    D2_m = foundation_results.D2_req_m

    if inputs.width_mode not in ("constant", "linear_to_D2"):
        raise ValueError("width_mode must be 'constant' or 'linear_to_D2'.")
    if ASD_m2 <= 0.0 or A2_m2 <= 0.0:
        raise ValueError("Areas must be positive.")

    # --------------------------------------------------------------------------
    # Width convention
    # --------------------------------------------------------------------------
    if inputs.width_mode == "constant":
        w2_eq_m = wSD_m
    else:
        # First-pass 2D-equivalent end-width tied to engine-face diameter
        w2_eq_m = D2_m

    h2_eq_m = A2_m2 / w2_eq_m

    # --------------------------------------------------------------------------
    # Candidate lengths
    # --------------------------------------------------------------------------
    L_area_based_m = compute_area_based_length(
        ASD_m2=ASD_m2,
        A2_m2=A2_m2,
        equivalent_conical_angle_deg=inputs.equivalent_conical_angle_deg,
    )

    L_width_based_m = compute_linear_dimension_length(
        start_val=wSD_m,
        end_val=w2_eq_m,
        equivalent_angle_deg=inputs.equivalent_width_transition_deg,
    )

    L_height_based_m = compute_linear_dimension_length(
        start_val=hSD_m,
        end_val=h2_eq_m,
        equivalent_angle_deg=inputs.equivalent_height_transition_deg,
    )

    # Choose the controlling length
    candidates = {
        "area_based": L_area_based_m,
        "width_based": L_width_based_m,
        "height_based": L_height_based_m,
    }
    controlling_mode = max(candidates, key=candidates.get)
    L_SD_to_2_m = candidates[controlling_mode]

    # --------------------------------------------------------------------------
    # Coordinates
    # Local coordinates for this script
    # --------------------------------------------------------------------------
    x_SD_m = 0.0
    x_2_m = L_SD_to_2_m

    # --------------------------------------------------------------------------
    # Ratios / interpretation
    # --------------------------------------------------------------------------
    A2_over_ASD = A2_m2 / ASD_m2
    h2_over_hSD = h2_eq_m / hSD_m
    w2_over_wSD = w2_eq_m / wSD_m

    percent_height_change = 100.0 * (h2_eq_m - hSD_m) / hSD_m
    percent_width_change = 100.0 * (w2_eq_m - wSD_m) / wSD_m

    interpretation_area = interpret_change(A2_m2, ASD_m2, "area")
    interpretation_width = interpret_change(w2_eq_m, wSD_m, "width")
    interpretation_height = interpret_change(h2_eq_m, hSD_m, "height")

    # --------------------------------------------------------------------------
    # Warnings
    # --------------------------------------------------------------------------
    if L_SD_to_2_m <= 0.0:
        warnings.append("Computed SD->2 length is non-positive.")

    if inputs.width_mode == "constant":
        warnings.append(
            "Width is held constant, so the rectangular-to-circular/annular transition is not yet geometrically represented."
        )
    else:
        warnings.append(
            "Width is linearly transitioned to D2 as a first-pass proxy for the rectangular-to-circular/annular transition."
        )

    if A2_over_ASD < 1.0:
        warnings.append("Station-2 area is smaller than SD area; this implies contraction rather than net diffusion.")
    if inputs.equivalent_conical_angle_deg > 4.0:
        warnings.append("Equivalent conical angle is relatively steep for a conservative first-pass subsonic diffuser.")
    if inputs.equivalent_conical_angle_deg < 1.5:
        warnings.append("Equivalent conical angle is very shallow and may imply a long diffuser.")

    if controlling_mode != "area_based":
        warnings.append(
            f"SD->2 length is controlled by {controlling_mode.replace('_', ' ')}, not by area change alone."
        )

    if abs(percent_width_change) > 20.0:
        warnings.append("Width transition is large and will strongly influence final diffuser geometry.")
    if abs(percent_height_change) > 20.0:
        warnings.append("Height transition is large and will strongly influence final diffuser geometry.")

    return SDTo2Results(
        ASD_m2=ASD_m2,
        hSD_m=hSD_m,
        wSD_m=wSD_m,
        A2_m2=A2_m2,
        D2_m=D2_m,
        h2_eq_m=h2_eq_m,
        w2_eq_m=w2_eq_m,
        L_area_based_m=L_area_based_m,
        L_width_based_m=L_width_based_m,
        L_height_based_m=L_height_based_m,
        L_SD_to_2_m=L_SD_to_2_m,
        controlling_mode=controlling_mode,
        x_SD_m=x_SD_m,
        x_2_m=x_2_m,
        A2_over_ASD=A2_over_ASD,
        h2_over_hSD=h2_over_hSD,
        w2_over_wSD=w2_over_wSD,
        percent_height_change=percent_height_change,
        percent_width_change=percent_width_change,
        interpretation_area=interpretation_area,
        interpretation_width=interpretation_width,
        interpretation_height=interpretation_height,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def print_sd_to_2_report(results: SDTo2Results, inputs: SDTo2Inputs) -> None:
    print("=" * 100)
    print("SCRIPT 6 — SD TO 2 V2")
    print("=" * 100)

    print("\nPROCEDURE CHECKPOINT")
    print("-" * 100)
    print("Current known station: SD")
    print("This script solves   : station 2 link in area/height/width/length terms")
    print("This script stops at : station 2")
    print("This is the section where the rectangular-to-circular/annular transition belongs.")

    print("\nASSUMPTIONS")
    print("-" * 100)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED SD BASIS")
    print("-" * 100)
    print(f"SD area                                         : {results.ASD_m2:.6f} m^2")
    print(f"SD equivalent height                            : {results.hSD_m:.6f} m")
    print(f"SD width                                        : {results.wSD_m:.6f} m")

    print("\nIMPORTED STATION-2 BASIS")
    print("-" * 100)
    print(f"Station-2 required annulus area                 : {results.A2_m2:.6f} m^2")
    print(f"Station-2 required D2                           : {results.D2_m:.6f} m")

    print("\nUSER-SET DIFFUSER FACTORS")
    print("-" * 100)
    print(f"Equivalent conical angle                        : {inputs.equivalent_conical_angle_deg:.6f} deg")
    print(f"Equivalent width-transition angle               : {inputs.equivalent_width_transition_deg:.6f} deg")
    print(f"Equivalent height-transition angle              : {inputs.equivalent_height_transition_deg:.6f} deg")
    print(f"Width mode                                      : {inputs.width_mode}")

    print("\nEQUIVALENT STATION-2 TARGET")
    print("-" * 100)
    print(f"Equivalent station-2 width                      : {results.w2_eq_m:.6f} m")
    print(f"Equivalent station-2 height                     : {results.h2_eq_m:.6f} m")

    print("\nCANDIDATE LENGTHS")
    print("-" * 100)
    print(f"Area-based length                               : {results.L_area_based_m:.6f} m")
    print(f"Width-based length                              : {results.L_width_based_m:.6f} m")
    print(f"Height-based length                             : {results.L_height_based_m:.6f} m")
    print(f"Controlling mode                                : {results.controlling_mode}")

    print("\nSD TO 2 LENGTH")
    print("-" * 100)
    print(f"SD x                                            : {results.x_SD_m:.6f} m")
    print(f"Station-2 x                                     : {results.x_2_m:.6f} m")
    print(f"SD->2 chosen length                             : {results.L_SD_to_2_m:.6f} m")

    print("\nAREA / HEIGHT / WIDTH PROGRESSION")
    print("-" * 100)
    print(f"Area ratio, A2/ASD                              : {results.A2_over_ASD:.6f}")
    print(f"Height ratio, h2/hSD                            : {results.h2_over_hSD:.6f}")
    print(f"Width ratio, w2/wSD                             : {results.w2_over_wSD:.6f}")
    print(f"Percent height change                           : {results.percent_height_change:.6f} %")
    print(f"Percent width change                            : {results.percent_width_change:.6f} %")
    print(f"Area interpretation                             : {results.interpretation_area}")
    print(f"Width interpretation                            : {results.interpretation_width}")
    print(f"Height interpretation                           : {results.interpretation_height}")

    print("\nHANDOVER")
    print("-" * 100)
    print("Main-path stations are now defined through station 2 in first-pass area/width/height/length terms.")
    print("Next step is to turn the SD->2 section into explicit geometry and integrate bifurcation.")

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

    # Script 4 basis
    shock_results = compute_terminal_normal_shock(
        external_results=external_results,
        mdot_kg_s=foundation_inputs.mdot_capture_kg_s,
    )

    # Script 5 basis
    ns_th_sd_inputs = NSTHSDInputs()
    ns_th_sd_results = compute_ns_th_sd(
        shock_results=shock_results,
        inputs=ns_th_sd_inputs,
        mdot_kg_s=foundation_inputs.mdot_capture_kg_s,
    )

    # Script 6 inputs
    sd_to_2_inputs = SDTo2Inputs()

    sd_to_2_results = compute_sd_to_2(
        ns_th_sd_results=ns_th_sd_results,
        foundation_results=foundation_results,
        inputs=sd_to_2_inputs,
    )

    print_sd_to_2_report(sd_to_2_results, sd_to_2_inputs)


if __name__ == "__main__":
    main()