#!/usr/bin/env python3
"""
Station-0 / Station-2 consistency check

Purpose
-------
1. Compute freestream ambient conditions at station 0 from ISA altitude
2. Compute required intake capture area from actual mass flow
3. Compute station-2 static state from Tt2, Pt2, M2
4. Compute required station-2 annulus flow area from continuity
5. Compute reference station-2 geometry from Rt and Rh/Rt
6. Compute required station-2 geometry that matches the thermodynamics
7. Compare reference vs required station-2 geometry

Inputs
------
- altitude_ft
- M0
- mdot_capture_kg_s
- mdot_2_kg_s
- Tt2_K
- Pt2_kPa
- M2
- Rt2_ref_m
- rh_over_rt
- capture_area_margin_fraction

Assumptions
-----------
- ISA atmosphere from altitude
- Station 0 = ambient freestream static conditions
- Actual mass flow is used for capture sizing
- mdot_capture = mdot_2 unless bleed / bypass is later introduced
- Station-2 static state is derived from Tt2, Pt2, M2
- Station-2 required flow area is computed from continuity
- Station-2 geometry is defined from Rt and Rh/Rt
- Required geometry is derived assuming the same Rh/Rt as the reference geometry
- Constant gamma = 1.4 and R = 287.05 J/kg/K
- Upstream capture area and downstream station-2 required flow area are separate quantities and are reported separately.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class Inputs:
    # Station 0: ambient freestream
    altitude_ft: float = 60000.0
    M0: float = 1.8
    mdot_capture_kg_s: float = 159.19
    engines_total: int = 4

    # Station 2: thermodynamic state
    mdot_2_kg_s: float = 159.19
    Tt2_K: float = 357.049
    Pt2_kPa: float = 39.169
    M2: float = 0.751

    # Station 2: reference geometry
    Rt2_ref_m: float = 0.840332733
    rh_over_rt: float = 0.3

    # Practical capture sizing
    capture_area_margin_fraction: float = 0.08

    # Example 2D capture widths
    example_widths_m: List[float] = field(
        default_factory=lambda: [2.0, 2.2, 2.4, 2.5]
    )

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Single design-point preliminary sizing only.",
            "Station-0 freestream static conditions are derived from ISA altitude.",
            "Actual mass flow is used for capture sizing.",
            "All captured mass flow is assumed to reach station 2.",
            "Station-2 static state is derived from Tt2, Pt2, and M2.",
            "Station-2 required flow area is computed from continuity.",
            "Reference station-2 geometry is defined from Rt and Rh/Rt.",
            "Required station-2 geometry is derived using the same Rh/Rt as the reference geometry.",
            "Perfect-gas model with constant gamma = 1.4 and R = 287.05 J/kg/K.",
            "No shock, recovery, throat, diffuser, distortion, or off-design modelling is included.",
            "Upstream capture area and downstream station-2 required flow area are separate quantities and are reported separately.",
        ]
    )


@dataclass
class Results:
    # Station 0 ambient
    T0_K: Optional[float]
    p0_kPa: Optional[float]
    rho0_kg_m3: Optional[float]
    a0_m_s: Optional[float]
    V0_m_s: Optional[float]

    # Capture area
    A_capture_ideal_m2: Optional[float]
    A_capture_working_m2: Optional[float]
    D_capture_equiv_ideal_m: Optional[float]
    D_capture_equiv_working_m: Optional[float]
    capture_examples: List[tuple[float, Optional[float], Optional[float]]]

    # Station 2 thermo
    T2_K: Optional[float]
    p2_kPa: Optional[float]
    rho2_kg_m3: Optional[float]
    a2_m_s: Optional[float]
    V2_m_s: Optional[float]
    A2_required_m2: Optional[float]

    # Station 2 reference geometry
    Rh2_ref_m: Optional[float]
    D2_ref_m: Optional[float]
    A2_gross_ref_m2: Optional[float]
    A2_annulus_ref_m2: Optional[float]

    # Station 2 required geometry
    Rt2_req_m: Optional[float]
    Rh2_req_m: Optional[float]
    D2_req_m: Optional[float]
    A2_gross_req_m2: Optional[float]
    A2_annulus_req_m2: Optional[float]

    # Comparison
    annulus_area_mismatch_ref_m2: Optional[float]
    annulus_area_mismatch_ref_percent: Optional[float]
    Rt_mismatch_m: Optional[float]
    Rh_mismatch_m: Optional[float]
    D2_mismatch_m: Optional[float]
    capture_to_annulus_req_ideal: Optional[float]
    capture_to_annulus_req_working: Optional[float]

    warnings: List[str]


# Gas constants
GAMMA = 1.4
R = 287.05


def fmt(x: Optional[float], ndp: int = 4) -> str:
    return "N/A" if x is None else f"{x:.{ndp}f}"


def density_from_pT(p_kPa: float, T_K: float) -> float:
    return (p_kPa * 1000.0) / (R * T_K)


def speed_of_sound(T_K: float) -> float:
    return math.sqrt(GAMMA * R * T_K)


def T_static_from_total(Tt_K: float, M: float) -> float:
    return Tt_K / (1.0 + 0.5 * (GAMMA - 1.0) * M**2)


def p_static_from_total(Pt_kPa: float, M: float) -> float:
    factor = (1.0 + 0.5 * (GAMMA - 1.0) * M**2) ** (GAMMA / (GAMMA - 1.0))
    return Pt_kPa / factor


def isa_from_altitude_ft(altitude_ft: float) -> tuple[float, float]:
    """
    Returns ISA ambient static temperature [K] and pressure [kPa].
    Valid here up to 32 km, which covers 60,000 ft.
    """
    h = altitude_ft * 0.3048
    g0 = 9.80665

    if h <= 11000.0:
        T = 288.15 - 0.0065 * h
        p = 101325.0 * (T / 288.15) ** (g0 / (R * 0.0065))
        return T, p / 1000.0

    if h <= 20000.0:
        T = 216.65
        p11 = 22632.06
        p = p11 * math.exp(-g0 * (h - 11000.0) / (R * T))
        return T, p / 1000.0

    if h <= 32000.0:
        T20 = 216.65
        p20 = 5474.889
        lapse = 0.001
        T = T20 + lapse * (h - 20000.0)
        p = p20 * (T / T20) ** (-g0 / (R * lapse))
        return T, p / 1000.0

    raise ValueError("ISA model only implemented up to 32 km.")


def compute(inputs: Inputs) -> Results:
    warnings: List[str] = []

    # -------------------------
    # Basic input checks
    # -------------------------
    if inputs.altitude_ft < 0:
        warnings.append("Altitude is negative.")
    if inputs.M0 <= 0:
        warnings.append("M0 must be > 0.")
    if inputs.mdot_capture_kg_s <= 0:
        warnings.append("mdot_capture_kg_s must be > 0.")
    if inputs.mdot_2_kg_s <= 0:
        warnings.append("mdot_2_kg_s must be > 0.")
    if inputs.Tt2_K <= 0:
        warnings.append("Tt2_K must be > 0.")
    if inputs.Pt2_kPa <= 0:
        warnings.append("Pt2_kPa must be > 0.")
    if not (0 < inputs.M2 < 1):
        warnings.append("M2 must be in the subsonic range (0, 1).")
    if inputs.Rt2_ref_m <= 0:
        warnings.append("Rt2_ref_m must be > 0.")
    if not (0 <= inputs.rh_over_rt < 1):
        warnings.append("rh_over_rt must be in [0, 1).")

    # -------------------------
    # Station 0 ambient
    # -------------------------
    T0 = p0 = rho0 = a0 = V0 = None
    try:
        T0, p0 = isa_from_altitude_ft(inputs.altitude_ft)
        rho0 = density_from_pT(p0, T0)
        a0 = speed_of_sound(T0)
        V0 = inputs.M0 * a0
    except Exception as e:
        warnings.append(f"Station 0 failed: {e}")

    # -------------------------
    # Capture area
    # -------------------------
    A_capture_ideal = A_capture_working = None
    D_capture_ideal = D_capture_working = None
    capture_examples: List[tuple[float, Optional[float], Optional[float]]] = []

    if rho0 is not None and V0 is not None and inputs.mdot_capture_kg_s > 0:
        A_capture_ideal = inputs.mdot_capture_kg_s / (rho0 * V0)
        A_capture_working = A_capture_ideal * (1.0 + inputs.capture_area_margin_fraction)
        D_capture_ideal = math.sqrt(4.0 * A_capture_ideal / math.pi)
        D_capture_working = math.sqrt(4.0 * A_capture_working / math.pi)

        for w in inputs.example_widths_m:
            if w > 0:
                capture_examples.append((w, A_capture_ideal / w, A_capture_working / w))
            else:
                capture_examples.append((w, None, None))
                warnings.append(f"Skipped non-positive capture width: {w}")
    else:
        warnings.append("Capture area not computed due to invalid station-0 inputs.")

    # -------------------------
    # Station 2 thermodynamics
    # -------------------------
    T2 = p2 = rho2 = a2 = V2 = A2_required = None

    if inputs.Tt2_K > 0 and inputs.Pt2_kPa > 0 and 0 < inputs.M2 < 1 and inputs.mdot_2_kg_s > 0:
        try:
            T2 = T_static_from_total(inputs.Tt2_K, inputs.M2)
            p2 = p_static_from_total(inputs.Pt2_kPa, inputs.M2)
            rho2 = density_from_pT(p2, T2)
            a2 = speed_of_sound(T2)
            V2 = inputs.M2 * a2
            A2_required = inputs.mdot_2_kg_s / (rho2 * V2)
        except Exception as e:
            warnings.append(f"Station 2 thermo failed: {e}")
    else:
        warnings.append("Station 2 thermo inputs invalid.")

    # -------------------------
    # Station 2 reference geometry
    # -------------------------
    Rh2_ref = D2_ref = A2_gross_ref = A2_annulus_ref = None

    if inputs.Rt2_ref_m > 0 and 0 <= inputs.rh_over_rt < 1:
        try:
            Rh2_ref = inputs.rh_over_rt * inputs.Rt2_ref_m
            D2_ref = 2.0 * inputs.Rt2_ref_m
            A2_gross_ref = math.pi * inputs.Rt2_ref_m**2
            A2_annulus_ref = math.pi * (inputs.Rt2_ref_m**2 - Rh2_ref**2)
        except Exception as e:
            warnings.append(f"Reference geometry failed: {e}")
    else:
        warnings.append("Reference geometry inputs invalid.")

    # -------------------------
    # Station 2 required geometry
    # -------------------------
    Rt2_req = Rh2_req = D2_req = A2_gross_req = A2_annulus_req = None

    if A2_required is not None and A2_required > 0 and 0 <= inputs.rh_over_rt < 1:
        try:
            lam = inputs.rh_over_rt
            Rt2_req = math.sqrt(A2_required / (math.pi * (1.0 - lam**2)))
            Rh2_req = lam * Rt2_req
            D2_req = 2.0 * Rt2_req
            A2_gross_req = math.pi * Rt2_req**2
            A2_annulus_req = math.pi * (Rt2_req**2 - Rh2_req**2)
        except Exception as e:
            warnings.append(f"Required geometry failed: {e}")
    else:
        warnings.append("Required geometry could not be computed.")

    # -------------------------
    # Comparison
    # -------------------------
    annulus_area_mismatch_ref = annulus_area_mismatch_ref_percent = None
    Rt_mismatch = Rh_mismatch = D2_mismatch = None
    capture_to_annulus_req_ideal = capture_to_annulus_req_working = None

    if A2_required is not None and A2_annulus_ref is not None and A2_annulus_ref > 0:
        annulus_area_mismatch_ref = A2_required - A2_annulus_ref
        annulus_area_mismatch_ref_percent = 100.0 * annulus_area_mismatch_ref / A2_annulus_ref

    if Rt2_req is not None and inputs.Rt2_ref_m > 0:
        Rt_mismatch = Rt2_req - inputs.Rt2_ref_m
    if Rh2_req is not None and Rh2_ref is not None:
        Rh_mismatch = Rh2_req - Rh2_ref
    if D2_req is not None and D2_ref is not None:
        D2_mismatch = D2_req - D2_ref

    if A_capture_ideal is not None and A2_annulus_req is not None and A2_annulus_req > 0:
        capture_to_annulus_req_ideal = A_capture_ideal / A2_annulus_req
    if A_capture_working is not None and A2_annulus_req is not None and A2_annulus_req > 0:
        capture_to_annulus_req_working = A_capture_working / A2_annulus_req

    return Results(
        T0_K=T0,
        p0_kPa=p0,
        rho0_kg_m3=rho0,
        a0_m_s=a0,
        V0_m_s=V0,
        A_capture_ideal_m2=A_capture_ideal,
        A_capture_working_m2=A_capture_working,
        D_capture_equiv_ideal_m=D_capture_ideal,
        D_capture_equiv_working_m=D_capture_working,
        capture_examples=capture_examples,
        T2_K=T2,
        p2_kPa=p2,
        rho2_kg_m3=rho2,
        a2_m_s=a2,
        V2_m_s=V2,
        A2_required_m2=A2_required,
        Rh2_ref_m=Rh2_ref,
        D2_ref_m=D2_ref,
        A2_gross_ref_m2=A2_gross_ref,
        A2_annulus_ref_m2=A2_annulus_ref,
        Rt2_req_m=Rt2_req,
        Rh2_req_m=Rh2_req,
        D2_req_m=D2_req,
        A2_gross_req_m2=A2_gross_req,
        A2_annulus_req_m2=A2_annulus_req,
        annulus_area_mismatch_ref_m2=annulus_area_mismatch_ref,
        annulus_area_mismatch_ref_percent=annulus_area_mismatch_ref_percent,
        Rt_mismatch_m=Rt_mismatch,
        Rh_mismatch_m=Rh_mismatch,
        D2_mismatch_m=D2_mismatch,
        capture_to_annulus_req_ideal=capture_to_annulus_req_ideal,
        capture_to_annulus_req_working=capture_to_annulus_req_working,
        warnings=warnings,
    )


def print_report(inputs: Inputs, results: Results) -> None:
    print("=" * 92)
    print("STATION-0 / STATION-2 CONSISTENCY CHECK")
    print("=" * 92)

    print("\nSTATION 0 (AMBIENT FREESTREAM)")
    print("-" * 92)
    print(f"Altitude (ISA):                  {inputs.altitude_ft:.1f} ft")
    print(f"M0:                              {inputs.M0:.4f}")
    print(f"T0:                              {fmt(results.T0_K, 4)} K")
    print(f"p0:                              {fmt(results.p0_kPa, 4)} kPa")
    print(f"rho0:                            {fmt(results.rho0_kg_m3, 6)} kg/m^3")
    print(f"a0:                              {fmt(results.a0_m_s, 4)} m/s")
    print(f"V0:                              {fmt(results.V0_m_s, 4)} m/s")

    print("\nCAPTURE AREA")
    print("-" * 92)
    print(f"mdot_capture:                    {inputs.mdot_capture_kg_s:.4f} kg/s")
    print(f"Ideal capture area:              {fmt(results.A_capture_ideal_m2, 4)} m^2")
    print(f"Working capture area:            {fmt(results.A_capture_working_m2, 4)} m^2")
    print(f"Ideal equiv. diameter:           {fmt(results.D_capture_equiv_ideal_m, 4)} m")
    print(f"Working equiv. diameter:         {fmt(results.D_capture_equiv_working_m, 4)} m")

    print("\nSTATION 2 THERMO")
    print("-" * 92)
    print(f"mdot2:                           {inputs.mdot_2_kg_s:.4f} kg/s")
    print(f"Tt2:                             {inputs.Tt2_K:.4f} K")
    print(f"Pt2:                             {inputs.Pt2_kPa:.4f} kPa")
    print(f"M2:                              {inputs.M2:.4f}")
    print(f"T2:                              {fmt(results.T2_K, 4)} K")
    print(f"p2:                              {fmt(results.p2_kPa, 4)} kPa")
    print(f"rho2:                            {fmt(results.rho2_kg_m3, 6)} kg/m^3")
    print(f"V2:                              {fmt(results.V2_m_s, 4)} m/s")
    print(f"Required flow area:              {fmt(results.A2_required_m2, 6)} m^2")

    print("\nREFERENCE STATION 2 GEOMETRY")
    print("-" * 92)
    print(f"Rt_ref:                          {inputs.Rt2_ref_m:.6f} m")
    print(f"Rh/Rt:                           {inputs.rh_over_rt:.6f}")
    print(f"Rh_ref:                          {fmt(results.Rh2_ref_m, 6)} m")
    print(f"D2_ref:                          {fmt(results.D2_ref_m, 6)} m")
    print(f"Gross face area ref:             {fmt(results.A2_gross_ref_m2, 6)} m^2")
    print(f"Annulus area ref:                {fmt(results.A2_annulus_ref_m2, 6)} m^2")

    print("\nREQUIRED STATION 2 GEOMETRY (DESIGN BASIS)")
    print("-" * 92)
    print(f"Rt_req:                          {fmt(results.Rt2_req_m, 6)} m")
    print(f"Rh_req:                          {fmt(results.Rh2_req_m, 6)} m")
    print(f"D2_req:                          {fmt(results.D2_req_m, 6)} m")
    print(f"Gross face area req:             {fmt(results.A2_gross_req_m2, 6)} m^2")
    print(f"Annulus area req:                {fmt(results.A2_annulus_req_m2, 6)} m^2")

    print("\nREFERENCE VS REQUIRED COMPARISON")
    print("-" * 92)
    print(f"Required - reference annulus:    {fmt(results.annulus_area_mismatch_ref_m2, 6)} m^2")
    print(f"Annulus mismatch percent:        {fmt(results.annulus_area_mismatch_ref_percent, 4)} %")
    print(f"Rt_req - Rt_ref:                 {fmt(results.Rt_mismatch_m, 6)} m")
    print(f"Rh_req - Rh_ref:                 {fmt(results.Rh_mismatch_m, 6)} m")
    print(f"D2_req - D2_ref:                 {fmt(results.D2_mismatch_m, 6)} m")

    print("\nAREA RATIOS USING REQUIRED GEOMETRY")
    print("-" * 92)
    print(f"Capture / annulus req (ideal):   {fmt(results.capture_to_annulus_req_ideal, 4)}")
    print(f"Capture / annulus req (working): {fmt(results.capture_to_annulus_req_working, 4)}")

    print("\nEXAMPLE 2D CAPTURE OPENINGS")
    print("-" * 92)
    print(f"{'Width (m)':>12} {'Ideal h (m)':>18} {'Working h (m)':>18}")
    for w, hi, hw in results.capture_examples:
        print(f"{w:12.4f} {fmt(hi, 4):>18} {fmt(hw, 4):>18}")

    print("\nTOTAL AIRCRAFT CAPTURE AREA")
    print("-" * 92)
    total_ideal = None if results.A_capture_ideal_m2 is None else results.A_capture_ideal_m2 * inputs.engines_total
    total_working = None if results.A_capture_working_m2 is None else results.A_capture_working_m2 * inputs.engines_total
    print(f"Total ideal capture area:        {fmt(total_ideal, 4)} m^2")
    print(f"Total working capture area:      {fmt(total_working, 4)} m^2")

    print("\nDESIGN ENDPOINT SUMMARY")
    print("-" * 92)
    print(f"Upstream ideal capture area:     {fmt(results.A_capture_ideal_m2, 6)} m^2")
    print(f"Upstream working capture area:   {fmt(results.A_capture_working_m2, 6)} m^2")
    print(f"Downstream required S2 area:     {fmt(results.A2_required_m2, 6)} m^2")
    print(f"Downstream required Rt:          {fmt(results.Rt2_req_m, 6)} m")
    print(f"Downstream required Rh:          {fmt(results.Rh2_req_m, 6)} m")
    print(f"Downstream required D2:          {fmt(results.D2_req_m, 6)} m")
    print(f"Downstream required gross area:  {fmt(results.A2_gross_req_m2, 6)} m^2")
    print(f"Downstream required annulus:     {fmt(results.A2_annulus_req_m2, 6)} m^2")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 92)
        for i, w in enumerate(results.warnings, start=1):
            print(f"{i:2d}. {w}")

    print("=" * 92)


def default_inputs() -> Inputs:
    return Inputs(
    )

def main() -> None:
    inputs = default_inputs()
    results = compute(inputs)
    print_report(inputs, results)



if __name__ == "__main__":
    main()