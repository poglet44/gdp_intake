#!/usr/bin/env python3
"""
prelim_intake_sizing.py

First-pass preliminary sizing script for a supersonic aircraft inlet.

Purpose
-------
Given a cruise condition and per-engine mass flow demand, this script computes:
1. Freestream properties
2. Required ideal capture area
3. Working capture area with design margin
4. Equivalent axisymmetric capture diameter
5. Example 2D intake opening dimensions
6. Engine-face gross / annulus areas
7. Area ratios useful for preliminary intake design

Current default values are based on the user's SST case:
- M = 1.6
- h = 60,000 ft
- T_inf = 216.65 K
- p_inf = 7.172 kPa
- m_dot = 159.19 kg/s per engine
- 4 engines total
- engine face diameter = 1.516185702 m

Notes
-----
- This is a first-pass capture sizing tool, not a full intake design code.
- It does NOT yet model shocks, recovery, distortion, bleed, variable geometry,
  or off-design behaviour.
- Capture area is driven mainly by required mass flow and freestream conditions.
- If ambient T and P are known from your engine deck / design point, use those.
- If hub_tip_ratio is None, annulus area is not calculated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional


# -----------------------------
# Configuration data structure
# -----------------------------
@dataclass
class IntakeSizingInputs:
    # Flight condition
    mach: float = 1.6
    altitude_ft: float = 60000.0

    # Freestream state at design point
    # Use known values from engine deck / mission analysis if available
    t_inf_k: float = 216.65
    p_inf_kpa: float = 7.172

    # Gas properties
    gamma: float = 1.4
    r_j_per_kgk: float = 287.05

    # Engine / intake requirements
    engines_total: int = 4
    mdot_per_engine_kg_s: float = 159.19
    engine_face_diameter_m: float = 1.516185702

    # Optional engine-face hub-tip ratio
    # Example: 0.30 or 0.35 if known. Set to None if unknown.
    hub_tip_ratio: Optional[float] = None

    # Intake sizing margin for first-pass practical area
    # 0.08 means +8%
    capture_area_margin_fraction: float = 0.08

    # Example 2D capture widths to inspect
    example_2d_widths_m: List[float] = field(
        default_factory=lambda: [2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
    )


# -----------------------------
# Result data structure
# -----------------------------
@dataclass
class IntakeSizingResults:
    # Freestream
    rho_inf_kg_m3: float
    a_inf_m_s: float
    v_inf_m_s: float
    mass_flux_kg_m2_s: float

    # Capture sizing
    capture_area_ideal_m2: float
    capture_area_working_m2: float
    axisymmetric_capture_diameter_ideal_m: float
    axisymmetric_capture_diameter_working_m: float

    # Engine-face geometry
    engine_face_gross_area_m2: float
    engine_face_hub_diameter_m: Optional[float]
    engine_face_annulus_area_m2: Optional[float]

    # Ratios
    capture_to_face_gross_ideal: float
    capture_to_face_gross_working: float
    capture_to_face_annulus_ideal: Optional[float]
    capture_to_face_annulus_working: Optional[float]

    # 2D examples
    two_d_examples: List[tuple[float, float, float]]


# -----------------------------
# Core calculations
# -----------------------------
def circle_area(diameter_m: float) -> float:
    return math.pi * diameter_m**2 / 4.0


def equivalent_diameter_from_area(area_m2: float) -> float:
    return math.sqrt(4.0 * area_m2 / math.pi)


def compute_freestream_density(p_kpa: float, t_k: float, r: float) -> float:
    """rho = p / (R T), with p converted from kPa to Pa."""
    return (p_kpa * 1000.0) / (r * t_k)


def compute_speed_of_sound(gamma: float, r: float, t_k: float) -> float:
    return math.sqrt(gamma * r * t_k)


def compute_intake_sizing(inputs: IntakeSizingInputs) -> IntakeSizingResults:
    # Freestream
    rho_inf = compute_freestream_density(
        p_kpa=inputs.p_inf_kpa,
        t_k=inputs.t_inf_k,
        r=inputs.r_j_per_kgk,
    )
    a_inf = compute_speed_of_sound(
        gamma=inputs.gamma,
        r=inputs.r_j_per_kgk,
        t_k=inputs.t_inf_k,
    )
    v_inf = inputs.mach * a_inf
    mass_flux = rho_inf * v_inf

    # Capture area per engine
    capture_area_ideal = inputs.mdot_per_engine_kg_s / mass_flux
    capture_area_working = capture_area_ideal * (1.0 + inputs.capture_area_margin_fraction)

    # Axisymmetric equivalent capture diameter
    axisym_d_ideal = equivalent_diameter_from_area(capture_area_ideal)
    axisym_d_working = equivalent_diameter_from_area(capture_area_working)

    # Engine face gross area
    face_gross_area = circle_area(inputs.engine_face_diameter_m)

    # Engine-face hub / annulus area
    if inputs.hub_tip_ratio is not None:
        hub_diameter = inputs.hub_tip_ratio * inputs.engine_face_diameter_m
        hub_area = circle_area(hub_diameter)
        annulus_area = face_gross_area - hub_area
    else:
        hub_diameter = None
        annulus_area = None

    # Ratios
    capture_to_face_gross_ideal = capture_area_ideal / face_gross_area
    capture_to_face_gross_working = capture_area_working / face_gross_area

    if annulus_area is not None and annulus_area > 0.0:
        capture_to_face_annulus_ideal = capture_area_ideal / annulus_area
        capture_to_face_annulus_working = capture_area_working / annulus_area
    else:
        capture_to_face_annulus_ideal = None
        capture_to_face_annulus_working = None

    # Example 2D openings: (width, ideal height, working height)
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
        capture_area_ideal_m2=capture_area_ideal,
        capture_area_working_m2=capture_area_working,
        axisymmetric_capture_diameter_ideal_m=axisym_d_ideal,
        axisymmetric_capture_diameter_working_m=axisym_d_working,
        engine_face_gross_area_m2=face_gross_area,
        engine_face_hub_diameter_m=hub_diameter,
        engine_face_annulus_area_m2=annulus_area,
        capture_to_face_gross_ideal=capture_to_face_gross_ideal,
        capture_to_face_gross_working=capture_to_face_gross_working,
        capture_to_face_annulus_ideal=capture_to_face_annulus_ideal,
        capture_to_face_annulus_working=capture_to_face_annulus_working,
        two_d_examples=two_d_examples,
    )


# -----------------------------
# Output formatting
# -----------------------------
def fmt(value: Optional[float], ndp: int = 4) -> str:
    if value is None:
        return "N/A"
    return f"{value:.{ndp}f}"


def print_report(inputs: IntakeSizingInputs, results: IntakeSizingResults) -> None:
    print("=" * 72)
    print("PRELIMINARY SUPERSONIC INTAKE SIZING REPORT")
    print("=" * 72)

    print("\nINPUTS")
    print("-" * 72)
    print(f"Cruise Mach                        : {inputs.mach:.4f}")
    print(f"Cruise altitude                    : {inputs.altitude_ft:.1f} ft")
    print(f"Freestream temperature             : {inputs.t_inf_k:.4f} K")
    print(f"Freestream pressure                : {inputs.p_inf_kpa:.4f} kPa")
    print(f"Per-engine mass flow               : {inputs.mdot_per_engine_kg_s:.4f} kg/s")
    print(f"Total number of engines            : {inputs.engines_total:d}")
    print(f"Engine-face diameter               : {inputs.engine_face_diameter_m:.6f} m")
    print(f"Hub-tip ratio                      : {fmt(inputs.hub_tip_ratio, 4)}")
    print(f"Capture area margin                : {100.0 * inputs.capture_area_margin_fraction:.2f} %")

    print("\nFREESTREAM")
    print("-" * 72)
    print(f"Density, rho_inf                   : {results.rho_inf_kg_m3:.6f} kg/m^3")
    print(f"Speed of sound, a_inf              : {results.a_inf_m_s:.4f} m/s")
    print(f"Flight speed, V_inf                : {results.v_inf_m_s:.4f} m/s")
    print(f"Mass flux, rho*V                   : {results.mass_flux_kg_m2_s:.4f} kg/(m^2 s)")

    print("\nCAPTURE AREA PER ENGINE")
    print("-" * 72)
    print(f"Ideal capture area                 : {results.capture_area_ideal_m2:.4f} m^2")
    print(f"Working capture area               : {results.capture_area_working_m2:.4f} m^2")

    print("\nAXISYMMETRIC EQUIVALENT CAPTURE SIZE PER ENGINE")
    print("-" * 72)
    print(f"Ideal equivalent diameter          : {results.axisymmetric_capture_diameter_ideal_m:.4f} m")
    print(f"Working equivalent diameter        : {results.axisymmetric_capture_diameter_working_m:.4f} m")

    print("\nENGINE-FACE AREA")
    print("-" * 72)
    print(f"Gross face area                    : {results.engine_face_gross_area_m2:.4f} m^2")
    print(f"Hub diameter                       : {fmt(results.engine_face_hub_diameter_m, 4)} m")
    print(f"Annulus area                       : {fmt(results.engine_face_annulus_area_m2, 4)} m^2")

    print("\nAREA RATIOS")
    print("-" * 72)
    print(f"Capture / gross face area (ideal)  : {results.capture_to_face_gross_ideal:.4f}")
    print(f"Capture / gross face area (working): {results.capture_to_face_gross_working:.4f}")
    print(f"Capture / annulus area (ideal)     : {fmt(results.capture_to_face_annulus_ideal, 4)}")
    print(f"Capture / annulus area (working)   : {fmt(results.capture_to_face_annulus_working, 4)}")

    print("\nEXAMPLE 2D OPENINGS PER ENGINE")
    print("-" * 72)
    print(f"{'Width (m)':>12} {'Ideal height (m)':>20} {'Working height (m)':>22}")
    for width, ideal_h, working_h in results.two_d_examples:
        print(f"{width:12.4f} {ideal_h:20.4f} {working_h:22.4f}")

    print("\nWHOLE AIRCRAFT TOTAL CAPTURE AREA")
    print("-" * 72)
    total_ideal = results.capture_area_ideal_m2 * inputs.engines_total
    total_working = results.capture_area_working_m2 * inputs.engines_total
    print(f"Total ideal capture area           : {total_ideal:.4f} m^2")
    print(f"Total working capture area         : {total_working:.4f} m^2")

    print("\nNOTES")
    print("-" * 72)
    print("* Capture area is primarily driven by required mass flow and freestream state.")
    print("* Therefore, for the same engine demand and cruise condition,")
    print("  2D bifurcated and axisymmetric inlets will need approximately")
    print("  the same capture area, but realised with different shapes.")
    print("* This script is only the first-pass sizing step.")
    print("* Next stage should include engine-face annulus definition,")
    print("  contraction ratios, shock system, pressure recovery,")
    print("  and off-design operability.")
    print("=" * 72)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    inputs = IntakeSizingInputs(
        mach=1.6,
        altitude_ft=60000.0,
        t_inf_k=216.65,
        p_inf_kpa=7.172,
        gamma=1.4,
        r_j_per_kgk=287.05,
        engines_total=4,
        mdot_per_engine_kg_s=159.19,
        engine_face_diameter_m=1.516185702,
        hub_tip_ratio=None,   # Change to 0.30 or 0.35 if known
        capture_area_margin_fraction=0.08,
        example_2d_widths_m=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5],
    )

    results = compute_intake_sizing(inputs)
    print_report(inputs, results)


if __name__ == "__main__":
    main()