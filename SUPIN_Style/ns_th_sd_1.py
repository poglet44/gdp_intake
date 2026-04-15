#!/usr/bin/env python3
"""
ns_th_sd_v1.py

Purpose
-------
Script 5 in the intake design chain.

This script solves the post-shock throat-section station structure:

    NS -> TH -> SD

It computes:
- TH state
- TH required area and equivalent height
- SD area and equivalent height
- simple section lengths NS->TH and TH->SD

It does NOT compute:
- detailed wall-line geometry
- bifurcation
- rectangular-to-circular transition
- SD->2 final diffuser geometry

Main path
---------
0 -> EX -> 1/NS -> TH -> SD -> 2
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


# --------------------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------------------

@dataclass
class NSTHSDInputs:
    """
    New design factors for NS -> TH -> SD.
    """
    MTH_target: float = 0.8
    ASD_over_ATH: float = 1.08
    theta_NS_to_TH_deg: float = 4.0
    theta_TH_to_SD_deg: float = 3.0

    assumptions: List[str] = field(
        default_factory=lambda: [
            "This script solves NS -> TH -> SD in state/area/length terms only.",
            "Width is held constant in this local 2D-equivalent model.",
            "TH state is computed from MTH_target using Tt_TH = Tt_NS and pt_TH = pt_NS in v1.",
            "No additional duct losses are included between NS and TH in this v1 model.",
            "SD area is set from the chosen area ratio ASD_over_ATH.",
            "Section lengths are computed using equivalent section-angle closure rules, not detailed wall-shape geometry.",
            "This script stops at SD and does not include SD->2, bifurcation, or rectangular-to-circular transition.",
        ]
    )


@dataclass
class NSTHSDResults:
    # Imported NS basis
    width_m: float
    MNS: float
    TNS_K: float
    pNS_kPa: float
    ptNS_kPa: float
    rhoNS_kg_m3: float
    aNS_m_s: float
    VNS_m_s: float
    ANS_required_m2: float
    hNS_m: float

    # TH state
    MTH: float
    TTH_K: float
    pTH_kPa: float
    ptTH_kPa: float
    rhoTH_kg_m3: float
    aTH_m_s: float
    VTH_m_s: float
    ATH_required_m2: float
    hTH_m: float

    # SD geometry/state placeholder
    ASD_required_m2: float
    hSD_m: float

    # Lengths
    L_NS_to_TH_m: float
    L_TH_to_SD_m: float
    L_NS_to_SD_m: float

    # Ratios / changes
    ATH_over_ANS: float
    ASD_over_ATH: float
    hTH_over_hNS: float
    hSD_over_hTH: float
    percent_height_change_NS_to_TH: float
    percent_height_change_TH_to_SD: float

    interpretation_NS_to_TH: str
    interpretation_TH_to_SD: str

    warnings: List[str]


# --------------------------------------------------------------------------------------
# Gas utilities
# --------------------------------------------------------------------------------------

def total_temperature_from_static(T_K: float, M: float, gamma: float) -> float:
    return T_K * (1.0 + 0.5 * (gamma - 1.0) * M * M)


def static_temperature_from_total(Tt_K: float, M: float, gamma: float) -> float:
    return Tt_K / (1.0 + 0.5 * (gamma - 1.0) * M * M)


def static_pressure_from_total(pt_kPa: float, M: float, gamma: float) -> float:
    return pt_kPa / (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------

def interpret_height_change(h2: float, h1: float) -> str:
    frac = (h2 - h1) / h1
    if frac > 0.15:
        return "aggressive expansion"
    if frac > 0.03:
        return "moderate expansion"
    if frac >= -0.03:
        return "near-constant height"
    if frac >= -0.15:
        return "moderate contraction"
    return "aggressive contraction"


def equivalent_length_from_height_change(
    h_start: float,
    h_end: float,
    theta_deg: float,
) -> float:
    """
    Convert height change into a first-pass section length using an equivalent
    section-angle closure rule.
    """
    if theta_deg <= 0.0:
        raise ValueError("Equivalent section angle must be > 0 deg.")

    delta_h = abs(h_end - h_start)
    if delta_h == 0.0:
        return 0.0

    return delta_h / tan(theta_deg * pi / 180.0)


# --------------------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------------------

def compute_ns_th_sd(
    shock_results,
    inputs: NSTHSDInputs,
    mdot_kg_s: float,
) -> NSTHSDResults:
    """
    Solve NS -> TH -> SD in state/area/length terms.
    """
    warnings: List[str] = []

    gamma = 1.4
    R = 287.05

    # --------------------------------------------------------------------------
    # Imported NS basis
    # --------------------------------------------------------------------------
    width = shock_results.capture_width_m
    MNS = shock_results.MNS
    TNS_K = shock_results.TNS_K
    pNS_kPa = shock_results.pNS_kPa
    ptNS_kPa = shock_results.ptNS_kPa
    rhoNS_kg_m3 = shock_results.rhoNS_kg_m3
    aNS_m_s = shock_results.aNS_m_s
    VNS_m_s = shock_results.VNS_m_s
    ANS_required_m2 = shock_results.ANS_required_m2
    hNS_m = shock_results.hNS_m

    # --------------------------------------------------------------------------
    # Input checks
    # --------------------------------------------------------------------------
    if not (MNS < inputs.MTH_target < 1.0):
        raise ValueError("MTH_target must satisfy MNS < MTH_target < 1.0.")
    if inputs.ASD_over_ATH <= 0.0:
        raise ValueError("ASD_over_ATH must be > 0.")
    if inputs.theta_NS_to_TH_deg <= 0.0:
        raise ValueError("theta_NS_to_TH_deg must be > 0.")
    if inputs.theta_TH_to_SD_deg <= 0.0:
        raise ValueError("theta_TH_to_SD_deg must be > 0.")

    # --------------------------------------------------------------------------
    # TH state from MTH_target
    # --------------------------------------------------------------------------
    TtNS_K = total_temperature_from_static(TNS_K, MNS, gamma)

    MTH = inputs.MTH_target
    TTH_K = static_temperature_from_total(TtNS_K, MTH, gamma)
    pTH_kPa = static_pressure_from_total(ptNS_kPa, MTH, gamma)
    rhoTH_kg_m3 = (pTH_kPa * 1000.0) / (R * TTH_K)
    aTH_m_s = sqrt(gamma * R * TTH_K)
    VTH_m_s = MTH * aTH_m_s
    ptTH_kPa = ptNS_kPa

    ATH_required_m2 = mdot_kg_s / (rhoTH_kg_m3 * VTH_m_s)
    hTH_m = ATH_required_m2 / width

    # --------------------------------------------------------------------------
    # SD area/height from area ratio
    # --------------------------------------------------------------------------
    ASD_required_m2 = inputs.ASD_over_ATH * ATH_required_m2
    hSD_m = ASD_required_m2 / width

    # --------------------------------------------------------------------------
    # Section lengths from equivalent angle rules
    # --------------------------------------------------------------------------
    L_NS_to_TH_m = equivalent_length_from_height_change(
        h_start=hNS_m,
        h_end=hTH_m,
        theta_deg=inputs.theta_NS_to_TH_deg,
    )

    L_TH_to_SD_m = equivalent_length_from_height_change(
        h_start=hTH_m,
        h_end=hSD_m,
        theta_deg=inputs.theta_TH_to_SD_deg,
    )

    L_NS_to_SD_m = L_NS_to_TH_m + L_TH_to_SD_m

    # --------------------------------------------------------------------------
    # Ratios / interpretation
    # --------------------------------------------------------------------------
    ATH_over_ANS = ATH_required_m2 / ANS_required_m2
    ASD_over_ATH = ASD_required_m2 / ATH_required_m2
    hTH_over_hNS = hTH_m / hNS_m
    hSD_over_hTH = hSD_m / hTH_m

    percent_height_change_NS_to_TH = 100.0 * (hTH_m - hNS_m) / hNS_m
    percent_height_change_TH_to_SD = 100.0 * (hSD_m - hTH_m) / hTH_m

    interpretation_NS_to_TH = interpret_height_change(hTH_m, hNS_m)
    interpretation_TH_to_SD = interpret_height_change(hSD_m, hTH_m)

    # --------------------------------------------------------------------------
    # Warnings
    # --------------------------------------------------------------------------
    if hTH_m <= 0.0:
        warnings.append("Computed TH opening height is non-positive.")
    if hSD_m <= 0.0:
        warnings.append("Computed SD opening height is non-positive.")
    if ATH_over_ANS > 1.0:
        warnings.append("TH area is larger than NS area; this is unusual for a throat section.")
    if ASD_over_ATH < 1.0:
        warnings.append("SD area is smaller than TH area; this is unusual if SD is the start of the subsonic diffuser.")
    if interpretation_NS_to_TH == "aggressive contraction":
        warnings.append("NS->TH contraction is aggressive and should be reviewed critically.")
    if interpretation_TH_to_SD == "aggressive expansion":
        warnings.append("TH->SD expansion is aggressive and should be reviewed critically.")

    return NSTHSDResults(
        width_m=width,
        MNS=MNS,
        TNS_K=TNS_K,
        pNS_kPa=pNS_kPa,
        ptNS_kPa=ptNS_kPa,
        rhoNS_kg_m3=rhoNS_kg_m3,
        aNS_m_s=aNS_m_s,
        VNS_m_s=VNS_m_s,
        ANS_required_m2=ANS_required_m2,
        hNS_m=hNS_m,
        MTH=MTH,
        TTH_K=TTH_K,
        pTH_kPa=pTH_kPa,
        ptTH_kPa=ptTH_kPa,
        rhoTH_kg_m3=rhoTH_kg_m3,
        aTH_m_s=aTH_m_s,
        VTH_m_s=VTH_m_s,
        ATH_required_m2=ATH_required_m2,
        hTH_m=hTH_m,
        ASD_required_m2=ASD_required_m2,
        hSD_m=hSD_m,
        L_NS_to_TH_m=L_NS_to_TH_m,
        L_TH_to_SD_m=L_TH_to_SD_m,
        L_NS_to_SD_m=L_NS_to_SD_m,
        ATH_over_ANS=ATH_over_ANS,
        ASD_over_ATH=ASD_over_ATH,
        hTH_over_hNS=hTH_over_hNS,
        hSD_over_hTH=hSD_over_hTH,
        percent_height_change_NS_to_TH=percent_height_change_NS_to_TH,
        percent_height_change_TH_to_SD=percent_height_change_TH_to_SD,
        interpretation_NS_to_TH=interpretation_NS_to_TH,
        interpretation_TH_to_SD=interpretation_TH_to_SD,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def print_ns_th_sd_report(results: NSTHSDResults, inputs: NSTHSDInputs) -> None:
    print("=" * 100)
    print("SCRIPT 5 — NS TO TH TO SD V1")
    print("=" * 100)

    print("\nPROCEDURE CHECKPOINT")
    print("-" * 100)
    print("Current known station: NS")
    print("This script solves   : TH and SD in state/area/length terms")
    print("This script stops at : SD")

    print("\nASSUMPTIONS")
    print("-" * 100)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED NS BASIS")
    print("-" * 100)
    print(f"Width                                           : {results.width_m:.6f} m")
    print(f"NS Mach                                         : {results.MNS:.6f}")
    print(f"NS static temperature                           : {results.TNS_K:.6f} K")
    print(f"NS static pressure                              : {results.pNS_kPa:.6f} kPa")
    print(f"NS total pressure                               : {results.ptNS_kPa:.6f} kPa")
    print(f"NS density                                      : {results.rhoNS_kg_m3:.6f} kg/m^3")
    print(f"NS speed of sound                               : {results.aNS_m_s:.6f} m/s")
    print(f"NS velocity                                     : {results.VNS_m_s:.6f} m/s")
    print(f"NS required area                                : {results.ANS_required_m2:.6f} m^2")
    print(f"NS opening height                               : {results.hNS_m:.6f} m")

    print("\nUSER-SET THROAT/DIFFUSER FACTORS")
    print("-" * 100)
    print(f"Target TH Mach, MTH_target                      : {inputs.MTH_target:.6f}")
    print(f"Area ratio, ASD/ATH                             : {inputs.ASD_over_ATH:.6f}")
    print(f"Equivalent angle, NS->TH                        : {inputs.theta_NS_to_TH_deg:.6f} deg")
    print(f"Equivalent angle, TH->SD                        : {inputs.theta_TH_to_SD_deg:.6f} deg")

    print("\nSOLVED TH STATE")
    print("-" * 100)
    print(f"TH Mach                                         : {results.MTH:.6f}")
    print(f"TH static temperature                           : {results.TTH_K:.6f} K")
    print(f"TH static pressure                              : {results.pTH_kPa:.6f} kPa")
    print(f"TH total pressure                               : {results.ptTH_kPa:.6f} kPa")
    print(f"TH density                                      : {results.rhoTH_kg_m3:.6f} kg/m^3")
    print(f"TH speed of sound                               : {results.aTH_m_s:.6f} m/s")
    print(f"TH velocity                                     : {results.VTH_m_s:.6f} m/s")
    print(f"TH required area                                : {results.ATH_required_m2:.6f} m^2")
    print(f"TH opening height                               : {results.hTH_m:.6f} m")

    print("\nSOLVED SD GEOMETRY TARGET")
    print("-" * 100)
    print(f"SD required area                                : {results.ASD_required_m2:.6f} m^2")
    print(f"SD opening height                               : {results.hSD_m:.6f} m")

    print("\nSECTION LENGTHS")
    print("-" * 100)
    print(f"NS to TH length                                 : {results.L_NS_to_TH_m:.6f} m")
    print(f"TH to SD length                                 : {results.L_TH_to_SD_m:.6f} m")
    print(f"NS to SD total length                           : {results.L_NS_to_SD_m:.6f} m")

    print("\nAREA / HEIGHT PROGRESSION")
    print("-" * 100)
    print(f"Area ratio, ATH/ANS                             : {results.ATH_over_ANS:.6f}")
    print(f"Area ratio, ASD/ATH                             : {results.ASD_over_ATH:.6f}")
    print(f"Height ratio, hTH/hNS                           : {results.hTH_over_hNS:.6f}")
    print(f"Height ratio, hSD/hTH                           : {results.hSD_over_hTH:.6f}")
    print(f"Percent height change NS->TH                    : {results.percent_height_change_NS_to_TH:.6f} %")
    print(f"Percent height change TH->SD                    : {results.percent_height_change_TH_to_SD:.6f} %")
    print(f"Interpretation NS->TH                           : {results.interpretation_NS_to_TH}")
    print(f"Interpretation TH->SD                           : {results.interpretation_TH_to_SD}")

    print("\nHANDOVER")
    print("-" * 100)
    print("Output stations solved: TH and SD")
    print("Next main script should start from SD and solve SD->2.")
    print("This script does not yet include detailed wall geometry, bifurcation, or rectangular-to-circular transition.")

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

    # Script 5 inputs
    ns_th_sd_inputs = NSTHSDInputs()

    ns_th_sd_results = compute_ns_th_sd(
        shock_results=shock_results,
        inputs=ns_th_sd_inputs,
        mdot_kg_s=foundation_inputs.mdot_capture_kg_s,
    )

    print_ns_th_sd_report(ns_th_sd_results, ns_th_sd_inputs)


if __name__ == "__main__":
    main()