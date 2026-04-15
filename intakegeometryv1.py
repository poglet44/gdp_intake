#!/usr/bin/env python3
"""
intake_geometry_v1.py

Purpose
-------
Complete Step 5: first-pass area-based intake geometry linking.

This script takes the two fixed endpoints from the foundation script:
- upstream capture requirement
- downstream required station-2 geometry

and produces a first-pass 2D bifurcated intake geometry summary.

What this script calculates
---------------------------
1. Capture-plane dimensions from chosen capture width
2. Symmetric branch flow-area split
3. First-pass area ratios between capture plane and station 2

Assumptions
-----------
1. The foundation script outputs are treated as fixed design-basis inputs.
2. Capture plane is represented as a simple rectangle:
       A_capture = width * height
3. Bifurcation is symmetric:
       each branch carries half the total flow area
4. Branch results are flow-area split only.
5. This script does NOT include:
   - splitter thickness
   - wall thickness
   - branch curvature
   - final bifurcation envelope
   - shocks / recovery / throat / diffuser design
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from intakesizingv2 import default_inputs, compute


@dataclass
class IntakeGeometryInputs:
    """
    Geometry choices introduced in Step 5.

    These are the only new choices allowed in this script.
    """
    capture_width_m: float = 1.85
    bifurcation_split_count: int = 2
    use_working_capture_area: bool = True

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Foundation-script outputs are treated as fixed design-basis inputs.",
            "Capture plane is represented as a simple width x height rectangle.",
            "Bifurcation is symmetric, so each branch carries half the total flow area.",
            "Branch results are flow-area split only, not final bifurcation geometry.",
            "No splitter thickness, wall thickness, curvature, shock system, throat, diffuser, or distortion modelling is included in this script.",
        ]
    )


@dataclass
class IntakeGeometryResults:
    # Imported fixed design basis
    A_capture_used_m2: float
    A2_annulus_req_m2: float
    A2_gross_req_m2: float
    Rt2_req_m: float
    Rh2_req_m: float
    D2_req_m: float

    # New Step 5 geometry results
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
    Compute Step 5 first-pass intake geometry.

    Method
    ------
    - import fixed upstream and downstream endpoint values
    - realise capture area using chosen width
    - split downstream required areas symmetrically into two branches
    - calculate simple area ratios
    """
    warnings: List[str] = []

    # -------------------------------------------------
    # Imported fixed design basis
    # -------------------------------------------------
    if geometry_inputs.use_working_capture_area:
        A_capture_used = foundation_results.A_capture_working_m2
    else:
        A_capture_used = foundation_results.A_capture_ideal_m2

    A2_annulus_req = foundation_results.A2_annulus_req_m2
    A2_gross_req = foundation_results.A2_gross_req_m2
    Rt2_req = foundation_results.Rt2_req_m
    Rh2_req = foundation_results.Rh2_req_m
    D2_req = foundation_results.D2_req_m

    # -------------------------------------------------
    # Checks
    # -------------------------------------------------
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

    # -------------------------------------------------
    # Capture-plane geometry
    # Assumption:
    # chosen capture width defines the rectangular opening height
    # -------------------------------------------------
    capture_height = A_capture_used / geometry_inputs.capture_width_m

    # -------------------------------------------------
    # First-pass symmetric bifurcation split
    # Assumption:
    # total downstream required area is split equally
    # -------------------------------------------------
    branch_annulus_area = A2_annulus_req / geometry_inputs.bifurcation_split_count
    branch_gross_area = A2_gross_req / geometry_inputs.bifurcation_split_count

    # -------------------------------------------------
    # Area ratios
    # -------------------------------------------------
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


def print_geometry_report(results: IntakeGeometryResults, geometry_inputs: IntakeGeometryInputs) -> None:
    """
    Print Step 5 report only.
    """
    print("=" * 92)
    print("STEP 5 — FIRST-PASS AREA-BASED INTAKE GEOMETRY")
    print("=" * 92)

    print("\nPURPOSE")
    print("-" * 92)
    print("Link the upstream capture requirement to the downstream required station-2 geometry.")
    print("This is a first-pass area-based geometry step only.")

    print("\nASSUMPTIONS")
    print("-" * 92)
    for i, assumption in enumerate(geometry_inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED FIXED DESIGN BASIS")
    print("-" * 92)
    print(f"Capture area used                              : {results.A_capture_used_m2:.6f} m^2")
    print(f"Required station-2 annulus area               : {results.A2_annulus_req_m2:.6f} m^2")
    print(f"Required station-2 gross area                 : {results.A2_gross_req_m2:.6f} m^2")
    print(f"Required station-2 Rt                         : {results.Rt2_req_m:.6f} m")
    print(f"Required station-2 Rh                         : {results.Rh2_req_m:.6f} m")
    print(f"Required station-2 D2                         : {results.D2_req_m:.6f} m")

    print("\nCAPTURE-PLANE GEOMETRY")
    print("-" * 92)
    print(f"Chosen capture width                          : {results.capture_width_m:.6f} m")
    print(f"Calculated capture height                     : {results.capture_height_m:.6f} m")

    print("\nFIRST-PASS BIFURCATION FLOW-AREA SPLIT")
    print("-" * 92)
    print(f"Symmetric split count                         : {geometry_inputs.bifurcation_split_count:d}")
    print(f"Branch annulus flow area per side            : {results.branch_annulus_area_m2:.6f} m^2")
    print(f"Branch gross area per side                   : {results.branch_gross_area_m2:.6f} m^2")

    print("\nAREA RATIOS")
    print("-" * 92)
    print(f"Total capture / required annulus              : {results.total_capture_to_annulus_ratio:.6f}")
    print(f"Total capture / required gross                : {results.total_capture_to_gross_ratio:.6f}")
    print(f"Per-branch capture / per-branch annulus       : {results.branch_capture_to_branch_annulus_ratio:.6f}")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 92)
        for i, warning in enumerate(results.warnings, start=1):
            print(f"{i:>2d}. {warning}")

    print("\nLIMITATION OF THIS STEP")
    print("-" * 92)
    print("This step does not define final bifurcation geometry.")
    print("It only defines first-pass flow-area based geometry quantities.")

    print("=" * 92)


def main() -> None:
    # -------------------------------------------------
    # Import fixed foundation design basis
    # -------------------------------------------------
    foundation_inputs = default_inputs()
    foundation_results = compute(foundation_inputs)

    # -------------------------------------------------
    # New geometry-only choices for Step 5
    # -------------------------------------------------
    geometry_inputs = IntakeGeometryInputs()

    geometry_results = compute_intake_geometry(foundation_results, geometry_inputs)
    print_geometry_report(geometry_results, geometry_inputs)


if __name__ == "__main__":
    main()