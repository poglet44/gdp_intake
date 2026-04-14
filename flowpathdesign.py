from intakesizingv2 import default_inputs, compute

def main() -> None:
    base_inputs = default_inputs()
    base_results = compute(base_inputs)

    A_capture = base_results.A_capture_working_m2
    D2_req = base_results.D2_req_m
    Rt2_req = base_results.Rt2_req_m
    Rh2_req = base_results.Rh2_req_m
    A2_req = base_results.A2_annulus_req_m2

    print(f"A_capture_working = {A_capture:.6f} m^2")
    print(f"D2_required       = {D2_req:.6f} m")
    print(f"Rt2_required      = {Rt2_req:.6f} m")
    print(f"Rh2_required      = {Rh2_req:.6f} m")
    print(f"A2_annulus_req    = {A2_req:.6f} m^2")

if __name__ == "__main__":
    main()