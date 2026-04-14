#!/usr/bin/env python3
"""
intake_geometry_v1.py

Purpose
-------
First-pass 2D bifurcated intake geometry sizing.

This script does NOT recalculate or redefine:
- station 0 ambient conditions
- capture area basis
- station 2 thermodynamic state
- required station-2 geometry

Those are imported from the foundation script and treated as fixed.

This script only adds:
- capture-plane geometry choice
- first-pass symmetric bifurcation split
- simple area-ratio geometry outputs
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

# Foundation script must expose:
# - default_inputs()
# - compute(inputs)
from intakesizingv2 import default_inputs, compute


@dataclass
class IntakeGeometryInputs:
    """
    New geometry-only inputs.

    These are the only values this script is allowed to choose.
    """
    capture_width_m: float = 1.85
    bifurcation_split_count: int = 2
    use_working_capture_area: bool = True


@dataclass
class IntakeGeometryResults:
    # Imported fixed design basis
    A_capture_used_m2: float
    A2_annulus_req_m2: float
    A2_gross_req_m2: float
    Rt2_req_m: float
    Rh2_req_m: float
    D2_req_m: float

    # New geometry results
    capture_width_m: float
    capture_height_m: float
    branch_annulus_area_m2: float
    branch_gross_area_m2: float

    # Ratios
    total_capture_to_annulus_ratio: float
    total_capture_to_gross_ratio: float
    branch_capture_to_branch_annulus_ratio: float

    warnings: List[str]


def compute_intake_geometry(
    foundation_results,
    geometry_inputs: IntakeGeometryInputs,
) -> IntakeGeometryResults:
    """
    Compute first-pass intake geometry using only imported foundation outputs.

    Assumptions
    -----------
    - Imported foundation results are treated as fixed design-basis inputs.
    - Capture plane is represented as a simple width x height rectangle.
    - Bifurcation is treated as a symmetric split.
    - Each branch carries equal area and equal flow at this stage.
    """
    warnings: List[str] = []

    # -------------------------
    # Fixed imported design basis
    # -------------------------
    if geometry_inputs.use_working_capture_area:
        A_capture_used = foundation_results.A_capture_working_m2
    else:
        A_capture_used = foundation_results.A_capture_ideal_m2

    A2_annulus_req = foundation_results.A2_annulus_req_m2
    A2_gross_req = foundation_results.A2_gross_req_m2
    Rt2_req = foundation_results.Rt2_req_m
    Rh2_req = foundation_results.Rh2_req_m
    D2_req = foundation_results.D2_req_m

    # -------------------------
    # Checks
    # -------------------------
    if A_capture_used is None or A_capture_used <= 0:
        raise ValueError("Imported capture area is invalid.")

    if A2_annulus_req is None or A2_annulus_req <= 0:
        raise ValueError("Imported required station-2 annulus area is invalid.")

    if A2_gross_req is None or A2_gross_req <= 0:
        raise ValueError("Imported required station-2 gross area is invalid.")

    if Rt2_req is None or Rh2_req is None or D2_req is None:
        raise ValueError("Imported required station-2 geometry is incomplete.")

    if geometry_inputs.capture_width_m <= 0:
        raise ValueError("capture_width_m must be > 0.")

    if geometry_inputs.bifurcation_split_count <= 0:
        raise ValueError("bifurcation_split_count must be > 0.")

    # -------------------------
    # Capture-plane geometry
    # -------------------------
    capture_height = A_capture_used / geometry_inputs.capture_width_m

    # -------------------------
    # First-pass symmetric bifurcation split
    # -------------------------
    branch_annulus_area = A2_annulus_req / geometry_inputs.bifurcation_split_count
    branch_gross_area = A2_gross_req / geometry_inputs.bifurcation_split_count

    # -------------------------
    # Ratios
    # -------------------------
    total_capture_to_annulus_ratio = A_capture_used / A2_annulus_req
    total_capture_to_gross_ratio = A_capture_used / A2_gross_req

    branch_capture_to_branch_annulus_ratio = (
        (A_capture_used / geometry_inputs.bifurcation_split_count) / branch_annulus_area
    )

    return IntakeGeometryResults(
        A_capture_used_m2=A_capture_used,
        A2_annulus_req_m2=A2_annulus_req,
        A2_gross_req_m2=A2_gross_req,
        Rt2_req_m=Rt2_req,
        Rh2_req_m=Rh2_req,
        D2_req_m=D2_req,
        capture_width_m=geometry_inputs.capture_width_m,
        capture_height_m=capture_height,
        branch_annulus_area_m2=branch_annulus_area,
        branch_gross_area_m2=branch_gross_area,
        total_capture_to_annulus_ratio=total_capture_to_annulus_ratio,
        total_capture_to_gross_ratio=total_capture_to_gross_ratio,
        branch_capture_to_branch_annulus_ratio=branch_capture_to_branch_annulus_ratio,
        warnings=warnings,
    )


def print_geometry_report(results: IntakeGeometryResults) -> None:
    """
    Print only the new geometry outputs.
    """
    print("=" * 88)
    print("FIRST-PASS 2D BIFURCATED INTAKE GEOMETRY")
    print("=" * 88)

    print("\nIMPORTED FIXED DESIGN BASIS")
    print("-" * 88)
    print(f"Capture area used                          : {results.A_capture_used_m2:.6f} m^2")
    print(f"Required station-2 annulus area           : {results.A2_annulus_req_m2:.6f} m^2")
    print(f"Required station-2 gross area             : {results.A2_gross_req_m2:.6f} m^2")
    print(f"Required station-2 Rt                     : {results.Rt2_req_m:.6f} m")
    print(f"Required station-2 Rh                     : {results.Rh2_req_m:.6f} m")
    print(f"Required station-2 D2                     : {results.D2_req_m:.6f} m")

    print("\nCAPTURE-PLANE GEOMETRY")
    print("-" * 88)
    print(f"Capture width                             : {results.capture_width_m:.6f} m")
    print(f"Capture height                            : {results.capture_height_m:.6f} m")

    print("\nFIRST-PASS BIFURCATION SPLIT")
    print("-" * 88)
    print(f"Branch annulus area per side             : {results.branch_annulus_area_m2:.6f} m^2")
    print(f"Branch gross area per side               : {results.branch_gross_area_m2:.6f} m^2")

    print("\nAREA RATIOS")
    print("-" * 88)
    print(f"Total capture / required annulus          : {results.total_capture_to_annulus_ratio:.6f}")
    print(f"Total capture / required gross            : {results.total_capture_to_gross_ratio:.6f}")
    print(f"Per-branch capture / per-branch annulus   : {results.branch_capture_to_branch_annulus_ratio:.6f}")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 88)
        for i, w in enumerate(results.warnings, start=1):
            print(f"{i:2d}. {w}")

    print("=" * 88)


def main() -> None:
    # -------------------------
    # Import fixed foundation design basis
    # -------------------------
    foundation_inputs = default_inputs()
    foundation_results = compute(foundation_inputs)

    # -------------------------
    # New geometry-only choices
    # -------------------------
    geometry_inputs = IntakeGeometryInputs()

    geometry_results = compute_intake_geometry(foundation_results, geometry_inputs)
    print_geometry_report(geometry_results)


if __name__ == "__main__":
    main()