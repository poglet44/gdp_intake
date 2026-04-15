#!/usr/bin/env python3
"""
intake_flowpath_progression_v1.py

Purpose
-------
First-pass intake flowpath progression between:
- upstream capture plane
- downstream required station-2 geometry

This script does NOT recalculate:
- station 0 ambient conditions
- capture area
- station 2 thermodynamics
- required station-2 geometry

Those are imported from the previous scripts and treated as fixed.

What this script does
---------------------
1. Defines pre-bifurcation and post-bifurcation lengths
2. Computes branch entry and exit flow areas
3. Computes first-pass rectangular branch dimensions
4. Computes simple contraction ratios and height-change rates

Assumptions
-----------
1. Imported endpoint values are fixed design-basis inputs.
2. Pre-bifurcation section is a single rectangular duct.
3. Bifurcation is symmetric.
4. No splitter thickness or wall thickness is included yet.
5. Branch width is chosen by the user.
6. Post-bifurcation branch width is held constant.
7. No shock, ramp, throat, diffuser, or off-design modelling is included.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from intakesizingv2 import default_inputs, compute
from intakegeometryv1 import IntakeGeometryInputs, compute_intake_geometry


@dataclass
class FlowpathInputs:
    """
    New geometry progression choices for this script only.
    """
    x_bifurcation_start_m: float = 4.0
    x_station2_m: float = 8.0
    branch_width_m: float = 1.2

    assumptions: List[str] = field(
        default_factory=lambda: [
            "Imported endpoint values are treated as fixed design-basis inputs.",
            "Pre-bifurcation section is represented as a single rectangular duct.",
            "Bifurcation is symmetric, so each branch carries half the total flow area.",
            "No splitter thickness or wall thickness is included yet.",
            "Branch width is user-chosen and held constant after the split.",
            "Branch heights are solved from flow area and chosen width.",
            "This script is geometry progression only; no ramps, shocks, throat, diffuser, or performance modelling are included.",
        ]
    )


@dataclass
class FlowpathResults:
    # Imported fixed basis
    capture_width_m: float
    capture_height_m: float
    A_capture_used_m2: float
    A2_annulus_req_m2: float
    A2_gross_req_m2: float
    Rt2_req_m: float
    Rh2_req_m: float
    D2_req_m: float

    # New progression results
    pre_bifurcation_length_m: float
    post_bifurcation_length_m: float

    branch_entry_area_m2: float
    branch_exit_area_m2: float

    branch_width_m: float
    branch_entry_height_m: float
    branch_exit_height_m: float

    branch_area_contraction_ratio: float
    branch_height_change_per_m: float

    warnings: List[str]


def compute_flowpath_progression(
    geometry_results,
    flowpath_inputs: FlowpathInputs,
) -> FlowpathResults:
    """
    Compute first-pass flowpath progression between capture plane and station 2.

    Method
    ------
    - import fixed Step 5 geometry basis
    - define pre- and post-bifurcation lengths
    - split capture and station-2 annulus areas symmetrically
    - compute branch heights from chosen branch width
    - compute contraction ratio and height gradient
    """
    warnings: List[str] = []

    # -------------------------------------------------
    # Imported fixed design basis
    # -------------------------------------------------
    capture_width = geometry_results.capture_width_m
    capture_height = geometry_results.capture_height_m
    A_capture_used = geometry_results.A_capture_used_m2

    A2_annulus_req = geometry_results.A2_annulus_req_m2
    A2_gross_req = geometry_results.A2_gross_req_m2
    Rt2_req = geometry_results.Rt2_req_m
    Rh2_req = geometry_results.Rh2_req_m
    D2_req = geometry_results.D2_req_m

    # -------------------------------------------------
    # Checks
    # -------------------------------------------------
    if flowpath_inputs.x_bifurcation_start_m <= 0:
        raise ValueError("x_bifurcation_start_m must be > 0.")

    if flowpath_inputs.x_station2_m <= flowpath_inputs.x_bifurcation_start_m:
        raise ValueError("x_station2_m must be greater than x_bifurcation_start_m.")

    if flowpath_inputs.branch_width_m <= 0:
        raise ValueError("branch_width_m must be > 0.")

    if A_capture_used is None or A_capture_used <= 0:
        raise ValueError("Imported capture area is invalid.")

    if A2_annulus_req is None or A2_annulus_req <= 0:
        raise ValueError("Imported required station-2 annulus area is invalid.")

    # -------------------------------------------------
    # Lengths
    # -------------------------------------------------
    L_pre = flowpath_inputs.x_bifurcation_start_m
    L_post = flowpath_inputs.x_station2_m - flowpath_inputs.x_bifurcation_start_m

    # -------------------------------------------------
    # Branch flow areas
    # Assumption:
    # symmetric bifurcation so each branch carries half the total flow area
    # -------------------------------------------------
    branch_entry_area = A_capture_used / 2.0
    branch_exit_area = A2_annulus_req / 2.0

    # -------------------------------------------------
    # Branch dimensions
    # Assumption:
    # branch width is fixed by the user
    # branch height changes to satisfy area change
    # -------------------------------------------------
    branch_width = flowpath_inputs.branch_width_m
    branch_entry_height = branch_entry_area / branch_width
    branch_exit_height = branch_exit_area / branch_width

    # -------------------------------------------------
    # Ratios and gradients
    # -------------------------------------------------
    branch_area_contraction_ratio = branch_entry_area / branch_exit_area
    branch_height_change_per_m = (branch_entry_height - branch_exit_height) / L_post

    return FlowpathResults(
        capture_width_m=capture_width,
        capture_height_m=capture_height,
        A_capture_used_m2=A_capture_used,
        A2_annulus_req_m2=A2_annulus_req,
        A2_gross_req_m2=A2_gross_req,
        Rt2_req_m=Rt2_req,
        Rh2_req_m=Rh2_req,
        D2_req_m=D2_req,
        pre_bifurcation_length_m=L_pre,
        post_bifurcation_length_m=L_post,
        branch_entry_area_m2=branch_entry_area,
        branch_exit_area_m2=branch_exit_area,
        branch_width_m=branch_width,
        branch_entry_height_m=branch_entry_height,
        branch_exit_height_m=branch_exit_height,
        branch_area_contraction_ratio=branch_area_contraction_ratio,
        branch_height_change_per_m=branch_height_change_per_m,
        warnings=warnings,
    )


def print_flowpath_report(results: FlowpathResults, flowpath_inputs: FlowpathInputs) -> None:
    print("=" * 92)
    print("STEP 6 — FIRST-PASS FLOWPATH GEOMETRY PROGRESSION")
    print("=" * 92)

    print("\nPURPOSE")
    print("-" * 92)
    print("Define the first-pass geometry progression between the capture plane and station 2.")
    print("This remains a geometry-only script.")

    print("\nASSUMPTIONS")
    print("-" * 92)
    for i, assumption in enumerate(flowpath_inputs.assumptions, start=1):
        print(f"{i:>2d}. {assumption}")

    print("\nIMPORTED FIXED DESIGN BASIS")
    print("-" * 92)
    print(f"Capture width                                  : {results.capture_width_m:.6f} m")
    print(f"Capture height                                 : {results.capture_height_m:.6f} m")
    print(f"Capture area used                              : {results.A_capture_used_m2:.6f} m^2")
    print(f"Required station-2 annulus area                : {results.A2_annulus_req_m2:.6f} m^2")
    print(f"Required station-2 gross area                  : {results.A2_gross_req_m2:.6f} m^2")
    print(f"Required station-2 Rt                          : {results.Rt2_req_m:.6f} m")
    print(f"Required station-2 Rh                          : {results.Rh2_req_m:.6f} m")
    print(f"Required station-2 D2                          : {results.D2_req_m:.6f} m")

    print("\nFLOWPATH LENGTHS")
    print("-" * 92)
    print(f"Pre-bifurcation length                         : {results.pre_bifurcation_length_m:.6f} m")
    print(f"Post-bifurcation length                        : {results.post_bifurcation_length_m:.6f} m")

    print("\nBRANCH FLOW AREAS")
    print("-" * 92)
    print(f"Branch entry flow area                         : {results.branch_entry_area_m2:.6f} m^2")
    print(f"Branch exit flow area                          : {results.branch_exit_area_m2:.6f} m^2")

    print("\nBRANCH RECTANGULAR DIMENSIONS")
    print("-" * 92)
    print(f"Chosen branch width                            : {results.branch_width_m:.6f} m")
    print(f"Branch entry height                            : {results.branch_entry_height_m:.6f} m")
    print(f"Branch exit height                             : {results.branch_exit_height_m:.6f} m")

    print("\nCONTRACTION METRICS")
    print("-" * 92)
    print(f"Branch area contraction ratio                  : {results.branch_area_contraction_ratio:.6f}")
    print(f"Branch height change per metre                 : {results.branch_height_change_per_m:.6f} m/m")

    if results.warnings:
        print("\nWARNINGS")
        print("-" * 92)
        for i, warning in enumerate(results.warnings, start=1):
            print(f"{i:>2d}. {warning}")

    print("\nLIMITATION OF THIS STEP")
    print("-" * 92)
    print("This step does not include splitter thickness, wall thickness, ramps, shocks,")
    print("throat sizing, diffuser performance, or final bifurcation geometry.")
    print("=" * 92)


def main() -> None:
    # -------------------------------------------------
    # Import fixed foundation design basis
    # -------------------------------------------------
    foundation_inputs = default_inputs()
    foundation_results = compute(foundation_inputs)

    # -------------------------------------------------
    # Import fixed Step 5 geometry basis
    # -------------------------------------------------
    geometry_inputs = IntakeGeometryInputs()
    geometry_results = compute_intake_geometry(foundation_results, geometry_inputs)

    # -------------------------------------------------
    # New geometry progression choices for this script
    # -------------------------------------------------
    flowpath_inputs = FlowpathInputs()

    flowpath_results = compute_flowpath_progression(geometry_results, flowpath_inputs)
    print_flowpath_report(flowpath_results, flowpath_inputs)


if __name__ == "__main__":
    main()