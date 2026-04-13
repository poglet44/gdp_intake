#!/usr/bin/env python3
"""
prelim_intake_sizing.py

First-pass preliminary sizing tool for a supersonic aircraft intake.

Current purpose
---------------
This version only performs capture-area-level sizing and basic geometry comparison.
It is intended for early concept work, not final intake validation.
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
#    No aircraft forebody, wing, or local installation effects are included yet.
#
# 3. Ideal capture area is based on continuity only:
#       A_capture = mdot / (rho_inf * V_inf)
#
# 4. No explicit modelling yet of:
#    - shocks
#    - total pressure recovery
#    - throat sizing
#    - diffuser losses
#    - bleed / bypass
#    - fan-face distortion
#    - 3D bifurcation effects
#    - off-design operability
#
# 5. A simple user-defined margin is applied to move from ideal capture area
#    to a practical "working" capture area.
#
# 6. For the same required mass flow and freestream condition, different inlet
#    families are assumed to need approximately the same capture area.
#    Therefore, at this stage, 2D bifurcated and axisymmetric concepts differ
#    mainly in how that area is realised geometrically.
#
# 7. Perfect-gas behaviour is assumed with constant gamma and R.
#
# 8. If hub-tip ratio is not supplied, annulus-based comparisons are omitted.
# =============================================================================


@dataclass
class IntakeSizingInputs:
    # -------------------------
    # Design point
    # -------------------------
    mach: float = 1.6
    altitude_ft: float = 60000.0

    # Freestream static conditions at design point
    # Replace with better values if you later change atmosphere source
    t_inf_k: float = 216.65
    p_inf_kpa: float = 7.172

    # Perfect-gas properties
    gamma: float = 1.4
    r_j_per_kgk: float = 287.05

    # -------------------------
    # Engine / intake inputs
    # -------------------------
    engines_total: int = 4
    mdot_per_engine_kg_s: float = 159.19
    engine_face_diameter_m: float = 1.516185702

    # Optional hub-tip ratio at engine face, e.g. 0.30 or 0.35
    hub_tip_ratio: Optional[float] = None

    # Simple practical sizing margin, e.g. 0.08 = +8%
    capture_area_margin_fraction: float = 0.08

    # Example 2D capture widths to inspect
    example_2d_widths_m: List[float] = field(
        default_factory=lambda: [2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
    )

    # Assumptions printed in the report
    assumptions: List[str] = field(
        default_factory=lambda: [
            "Single design-point preliminary sizing only.",
            "Uniform freestream assumed at the intake capture plane.",
            "Ideal capture area computed from continuity only: A = mdot / (rho * V).",
            "No shock, recovery, throat, diffuser, bleed, bypass, or distortion modelling in this version.",
            "A simple sizing margin is used to move from ideal to working capture area.",
            "2D bifurcated and axisymmetric concepts are assumed to require approximately the same capture area for the same mdot and freestream state.",
            "Perfect-gas model with constant gamma and gas constant R.",
            "Annulus-based comparisons are omitted if hub-tip ratio is not supplied.",
        ]
    )


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

    # 2D examples: (width, ideal height, working height)
    two_d_examples: List[tuple[float, float, float]]


def circle_area(diameter_m: float) -> float:
    """Return circular area from diameter."""
    return math.pi * diameter_m**2 / 4.0


def equivalent_diameter_from_area(area_m2: float) -> float:
    """Return equivalent circular diameter for a given area."""
    return math.sqrt(4.0 * area_m2 / math.pi)


def compute_freestream_density(p_kpa: float, t_k: float, r: float) -> float:
    """Compute density using rho = p / (R T), converting kPa to Pa."""
    return (p_kpa * 1000.0) / (r * t_k)


def compute_speed_of_sound(gamma: float, r: float, t_k: float) -> float:
    """Compute speed of sound for a perfect gas."""
    return math.sqrt(gamma * r * t_k)


def compute_intake_sizing(inputs: IntakeSizingInputs) -> IntakeSizingResults:
    """Run first-pass intake sizing calculations."""

    # Freestream properties
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

    # Equivalent axisymmetric capture diameter
    axisym_d_ideal = equivalent_diameter_from_area(capture_area_ideal)
    axisym_d_working = equivalent_diameter_from_area(capture_area_working)

    # Engine-face gross area
    face_gross_area = circle_area(inputs.engine_face_diameter_m)

    # Engine-face annulus area if hub-tip ratio is known
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

    # Example 2D openings
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


def fmt(value: Optional[float], ndp: int = 4) -> str:
    """Format optional float values for report printing."""
    if value is None:
        return "N/A"
    return f"{value:.{ndp}f}"


def print_report(inputs: IntakeSizingInputs, results: IntakeSizingResults) -> None:
    """Print a clear summary report to console."""

    print("=" * 80)
    print("PRELIMINARY SUPERSONIC INTAKE SIZING REPORT")
    print("=" * 80)

    print("\nPURPOSE")
    print("-" * 80)
    print("First-pass intake sizing at a single design point for early concept work.")

    print("\nINPUTS")
    print("-" * 80)
    print(f"Cruise Mach                                : {inputs.mach:.4f}")
    print(f"Cruise altitude                            : {inputs.altitude_ft:.1f} ft")
    print(f"Freestream temperature                     : {inputs.t_inf_k:.4f} K")
    print(f"Freestream pressure                        : {inputs.p_inf_kpa:.4f} kPa")
    print(f"Per-engine mass flow                       : {inputs.mdot_per_engine_kg_s:.4f} kg/s")
    print(f"Total number of engines                    : {inputs.engines_total:d}")
    print(f"Engine-face diameter                       : {inputs.engine_face_diameter_m:.6f} m")
    print(f"Hub-tip ratio                              : {fmt(inputs.hub_tip_ratio, 4)}")
    print(f"Capture area margin                        : {100.0 * inputs.capture_area_margin_fraction:.2f} %")

    print("\nASSUMPTIONS")
    print("-" * 80)
    for i, assumption in enumerate(inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nFREESTREAM")
    print("-" * 80)
    print(f"Density, rho_inf                           : {results.rho_inf_kg_m3:.6f} kg/m^3")
    print(f"Speed of sound, a_inf                      : {results.a_inf_m_s:.4f} m/s")
    print(f"Flight speed, V_inf                        : {results.v_inf_m_s:.4f} m/s")
    print(f"Mass flux, rho_inf * V_inf                 : {results.mass_flux_kg_m2_s:.4f} kg/(m^2 s)")

    print("\nCAPTURE AREA PER ENGINE")
    print("-" * 80)
    print(f"Ideal capture area                         : {results.capture_area_ideal_m2:.4f} m^2")
    print(f"Working capture area                       : {results.capture_area_working_m2:.4f} m^2")

    print("\nAXISYMMETRIC EQUIVALENT CAPTURE SIZE PER ENGINE")
    print("-" * 80)
    print(f"Ideal equivalent capture diameter          : {results.axisymmetric_capture_diameter_ideal_m:.4f} m")
    print(f"Working equivalent capture diameter        : {results.axisymmetric_capture_diameter_working_m:.4f} m")

    print("\nENGINE FACE AREA")
    print("-" * 80)
    print(f"Gross engine-face area                     : {results.engine_face_gross_area_m2:.4f} m^2")
    print(f"Hub diameter                               : {fmt(results.engine_face_hub_diameter_m, 4)} m")
    print(f"Annulus area                               : {fmt(results.engine_face_annulus_area_m2, 4)} m^2")

    print("\nAREA RATIOS")
    print("-" * 80)
    print(f"Capture / gross face area (ideal)          : {results.capture_to_face_gross_ideal:.4f}")
    print(f"Capture / gross face area (working)        : {results.capture_to_face_gross_working:.4f}")
    print(f"Capture / annulus area (ideal)             : {fmt(results.capture_to_face_annulus_ideal, 4)}")
    print(f"Capture / annulus area (working)           : {fmt(results.capture_to_face_annulus_working, 4)}")

    print("\nEXAMPLE 2D OPENINGS PER ENGINE")
    print("-" * 80)
    print(f"{'Width (m)':>12} {'Ideal height (m)':>20} {'Working height (m)':>22}")
    for width, ideal_h, working_h in results.two_d_examples:
        print(f"{width:12.4f} {ideal_h:20.4f} {working_h:22.4f}")

    print("\nWHOLE AIRCRAFT TOTAL CAPTURE AREA")
    print("-" * 80)
    total_ideal = results.capture_area_ideal_m2 * inputs.engines_total
    total_working = results.capture_area_working_m2 * inputs.engines_total
    print(f"Total ideal capture area                   : {total_ideal:.4f} m^2")
    print(f"Total working capture area                 : {total_working:.4f} m^2")

    print("\nLIMITATION OF THIS VERSION")
    print("-" * 80)
    print("This script only handles first-pass capture sizing.")
    print("It does NOT yet predict pressure recovery, throat sizing,")
    print("subsonic diffuser performance, distortion, or off-design behaviour.")

    print("=" * 80)


def main() -> None:
    inputs = IntakeSizingInputs(
        mach=1.6,
        altitude_ft=60000.0,
        t_inf_k=216.65,
        p_inf_kpa=7.172,
        gamma=1.4,
        r_j_per_kgk=287.05,
        engines_total=4,
        mdot_per_engine_kg_s=159.19,   # Mass flow per engine at cruise
        engine_face_diameter_m=1.516185702,
        hub_tip_ratio=0.45,  
        capture_area_margin_fraction=0.08,
        example_2d_widths_m=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5],
    )

    results = compute_intake_sizing(inputs)
    print_report(inputs, results)


if __name__ == "__main__":
    main()