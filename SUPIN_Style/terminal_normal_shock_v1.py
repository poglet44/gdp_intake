#!/usr/bin/env python3
"""
terminal_normal_shock_v1.py

Purpose
-------
Script 4 in the intake design chain.

This script takes the EX state from Script 2 and applies the terminal normal shock
to compute the first real post-shock internal starting state, labelled NS.

This script does NOT:
- define throat geometry
- define bifurcation geometry
- define subsonic diffuser geometry
- include bleed / bypass / viscous corrections

It only computes:
- EX state (echoed)
- NS state after terminal normal shock
- required NS area and height from continuity
- EX-to-NS comparisons

Design meaning
--------------
This is the first script that tells you what the internal duct actually receives
after the terminal shock.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt
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
class TerminalShockResults:
    # Imported EX basis
    capture_width_m: float
    MEX: float
    TEX_K: float
    pEX_kPa: float
    ptEX_kPa: float
    rhoEX_kg_m3: float
    aEX_m_s: float
    VEX_m_s: float
    AEX_required_m2: float
    hEX_m: float

    # NS state
    MNS: float
    TNS_K: float
    pNS_kPa: float
    ptNS_kPa: float
    rhoNS_kg_m3: float
    aNS_m_s: float
    VNS_m_s: float
    ANS_required_m2: float
    hNS_m: float

    # Comparisons
    pNS_over_pEX: float
    ptNS_over_ptEX: float
    TNS_over_TEX: float
    ANS_over_AEX: float
    hNS_over_hEX: float
    delta_h_m: float
    percent_height_change: float
    interpretation: str

    warnings: List[str]


# --------------------------------------------------------------------------------------
# Normal shock utilities
# --------------------------------------------------------------------------------------

def normal_shock_M2(M1: float, gamma: float) -> float:
    return sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * M1 * M1) /
        (gamma * M1 * M1 - 0.5 * (gamma - 1.0))
    )


def normal_shock_p2_over_p1(M1: float, gamma: float) -> float:
    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (M1 * M1 - 1.0)


def normal_shock_rho2_over_rho1(M1: float, gamma: float) -> float:
    return ((gamma + 1.0) * M1 * M1) / ((gamma - 1.0) * M1 * M1 + 2.0)


def normal_shock_T2_over_T1(M1: float, gamma: float) -> float:
    p_ratio = normal_shock_p2_over_p1(M1, gamma)
    rho_ratio = normal_shock_rho2_over_rho1(M1, gamma)
    return p_ratio / rho_ratio


def total_pressure_from_static(p_kPa: float, M: float, gamma: float) -> float:
    return p_kPa * (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


# --------------------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------------------

def compute_terminal_normal_shock(external_results, mdot_kg_s: float) -> TerminalShockResults:
    """
    Apply a terminal normal shock to the EX state from Script 2.
    """
    warnings: List[str] = []

    gamma = 1.4
    R = 287.05

    # Imported EX basis
    width = external_results.capture_width_m
    MEX = external_results.MEX_actual
    TEX_K = external_results.TEX_K
    pEX_kPa = external_results.pEX_kPa
    ptEX_kPa = external_results.ptEX_kPa
    rhoEX_kg_m3 = external_results.rhoEX_kg_m3
    aEX_m_s = external_results.aEX_m_s
    VEX_m_s = external_results.VEX_m_s
    AEX_required_m2 = external_results.AEX_required_m2
    hEX_m = external_results.hEX_m

    if MEX <= 1.0:
        raise ValueError("EX Mach must be > 1 for a terminal normal shock.")

    # Normal shock
    MNS = normal_shock_M2(MEX, gamma)
    pNS_over_pEX = normal_shock_p2_over_p1(MEX, gamma)
    rhoNS_over_rhoEX = normal_shock_rho2_over_rho1(MEX, gamma)
    TNS_over_TEX = normal_shock_T2_over_T1(MEX, gamma)

    pNS_kPa = pEX_kPa * pNS_over_pEX
    TNS_K = TEX_K * TNS_over_TEX
    rhoNS_kg_m3 = rhoEX_kg_m3 * rhoNS_over_rhoEX

    aNS_m_s = sqrt(gamma * R * TNS_K)
    VNS_m_s = MNS * aNS_m_s

    # Recompute downstream total pressure from downstream static state
    ptNS_kPa = total_pressure_from_static(pNS_kPa, MNS, gamma)
    ptNS_over_ptEX = ptNS_kPa / ptEX_kPa

    # Continuity
    ANS_required_m2 = mdot_kg_s / (rhoNS_kg_m3 * VNS_m_s)
    hNS_m = ANS_required_m2 / width

    ANS_over_AEX = ANS_required_m2 / AEX_required_m2
    hNS_over_hEX = hNS_m / hEX_m
    delta_h_m = hNS_m - hEX_m
    percent_height_change = 100.0 * delta_h_m / hEX_m

    # Interpretation
    if hNS_over_hEX < 1.05:
        interpretation = "mild post-shock area increase"
    elif hNS_over_hEX < 1.20:
        interpretation = "moderate post-shock area increase"
    else:
        interpretation = "large post-shock area increase"

    # Warnings
    if MNS >= 1.0:
        warnings.append("Post-shock Mach is not subsonic. Check the shock calculation.")
    if ptNS_over_ptEX >= 1.0:
        warnings.append("Total pressure did not decrease across the normal shock. Check the calculation.")
    if hNS_m <= 0.0:
        warnings.append("Computed NS height is non-positive.")
    if hNS_over_hEX > 1.25:
        warnings.append("Post-shock opening increase is large and may strongly affect internal geometry.")

    return TerminalShockResults(
        capture_width_m=width,
        MEX=MEX,
        TEX_K=TEX_K,
        pEX_kPa=pEX_kPa,
        ptEX_kPa=ptEX_kPa,
        rhoEX_kg_m3=rhoEX_kg_m3,
        aEX_m_s=aEX_m_s,
        VEX_m_s=VEX_m_s,
        AEX_required_m2=AEX_required_m2,
        hEX_m=hEX_m,
        MNS=MNS,
        TNS_K=TNS_K,
        pNS_kPa=pNS_kPa,
        ptNS_kPa=ptNS_kPa,
        rhoNS_kg_m3=rhoNS_kg_m3,
        aNS_m_s=aNS_m_s,
        VNS_m_s=VNS_m_s,
        ANS_required_m2=ANS_required_m2,
        hNS_m=hNS_m,
        pNS_over_pEX=pNS_over_pEX,
        ptNS_over_ptEX=ptNS_over_ptEX,
        TNS_over_TEX=TNS_over_TEX,
        ANS_over_AEX=ANS_over_AEX,
        hNS_over_hEX=hNS_over_hEX,
        delta_h_m=delta_h_m,
        percent_height_change=percent_height_change,
        interpretation=interpretation,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def print_terminal_shock_report(results: TerminalShockResults) -> None:
    print("=" * 100)
    print("SCRIPT 4 — TERMINAL NORMAL SHOCK / NS STATE V1")
    print("=" * 100)

    print("\nPURPOSE")
    print("-" * 100)
    print("Apply the terminal normal shock to the EX state from Script 2 to obtain the first")
    print("real post-shock internal starting state, NS.")

    print("\nIMPORTED EX BASIS")
    print("-" * 100)
    print(f"Width                                           : {results.capture_width_m:.6f} m")
    print(f"EX Mach                                         : {results.MEX:.6f}")
    print(f"EX static temperature                           : {results.TEX_K:.6f} K")
    print(f"EX static pressure                              : {results.pEX_kPa:.6f} kPa")
    print(f"EX total pressure                               : {results.ptEX_kPa:.6f} kPa")
    print(f"EX density                                      : {results.rhoEX_kg_m3:.6f} kg/m^3")
    print(f"EX speed of sound                               : {results.aEX_m_s:.6f} m/s")
    print(f"EX velocity                                     : {results.VEX_m_s:.6f} m/s")
    print(f"EX required area                                : {results.AEX_required_m2:.6f} m^2")
    print(f"EX opening height                               : {results.hEX_m:.6f} m")

    print("\nNS STATE (JUST DOWNSTREAM OF TERMINAL NORMAL SHOCK)")
    print("-" * 100)
    print(f"NS Mach                                         : {results.MNS:.6f}")
    print(f"NS static temperature                           : {results.TNS_K:.6f} K")
    print(f"NS static pressure                              : {results.pNS_kPa:.6f} kPa")
    print(f"NS total pressure                               : {results.ptNS_kPa:.6f} kPa")
    print(f"NS density                                      : {results.rhoNS_kg_m3:.6f} kg/m^3")
    print(f"NS speed of sound                               : {results.aNS_m_s:.6f} m/s")
    print(f"NS velocity                                     : {results.VNS_m_s:.6f} m/s")
    print(f"NS required area                                : {results.ANS_required_m2:.6f} m^2")
    print(f"NS opening height                               : {results.hNS_m:.6f} m")

    print("\nEX TO NS COMPARISON")
    print("-" * 100)
    print(f"Static pressure ratio, pNS/pEX                  : {results.pNS_over_pEX:.6f}")
    print(f"Total pressure ratio, ptNS/ptEX                 : {results.ptNS_over_ptEX:.6f}")
    print(f"Static temperature ratio, TNS/TEX               : {results.TNS_over_TEX:.6f}")
    print(f"Area ratio, ANS/AEX                             : {results.ANS_over_AEX:.6f}")
    print(f"Height ratio, hNS/hEX                           : {results.hNS_over_hEX:.6f}")
    print(f"Height change                                   : {results.delta_h_m:.6f} m")
    print(f"Percent height change                           : {results.percent_height_change:.6f} %")
    print(f"Interpretation                                  : {results.interpretation}")

    print("\nHOW TO TAKE THIS RESULT")
    print("-" * 100)
    print("This tells you how much the terminal normal shock changes the internal starting state.")
    print("The key design takeaway is the change from EX opening/area to NS opening/area.")

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

    # Script 4
    shock_results = compute_terminal_normal_shock(
        external_results=external_results,
        mdot_kg_s=foundation_inputs.mdot_capture_kg_s,
    )
    print_terminal_shock_report(shock_results)


if __name__ == "__main__":
    main()