"""Re-run the angle sweep with cascading warm starts from the
well-optimised alpha = 90 degree case outward in both directions."""

from __future__ import annotations

import math
import time
import numpy as np

from sofa_angle import optimise_sofa_alpha


ALPHAS_DEG = [60, 75, 90, 105, 120, 135, 150]


def run(out_path: str = "sweep_angles.npz") -> None:
    # Center anchor at 90 degrees, cascade outward
    results = {ad: None for ad in ALPHAS_DEG}
    t0 = time.time()

    print("=== 90 deg from quarter-arc seed ===", flush=True)
    results[90] = optimise_sofa_alpha(
        alpha=math.radians(90), M=10, n_thetas=80, n_rounds=3,
        maxiter=2500, seed=7)

    # Going up: 90 -> 105 -> 120 -> 135 -> 150
    prev = 90
    for ad in [105, 120, 135, 150]:
        print(f"=== {ad} deg, warm-start from {prev} ===", flush=True)
        results[ad] = optimise_sofa_alpha(
            alpha=math.radians(ad), M=10, n_thetas=80, n_rounds=4,
            maxiter=2500, seed=ad,
            warm_start=results[prev]["controls"])
        prev = ad

    # Going down: 90 -> 75 -> 60
    prev = 90
    for ad in [75, 60]:
        print(f"=== {ad} deg, warm-start from {prev} ===", flush=True)
        results[ad] = optimise_sofa_alpha(
            alpha=math.radians(ad), M=10, n_thetas=80, n_rounds=4,
            maxiter=2500, seed=ad,
            warm_start=results[prev]["controls"])
        prev = ad

    print(f"\ntotal time: {time.time() - t0:.0f}s", flush=True)
    print("\nfinal areas:")
    for ad in ALPHAS_DEG:
        print(f"  alpha={ad:3d}: area={results[ad]['area']:.5f}", flush=True)

    # Save bundle
    save_dict = {"alphas_deg": np.array(ALPHAS_DEG)}
    areas = []
    for i, ad in enumerate(ALPHAS_DEG):
        r = results[ad]
        areas.append(r["area"])
        save_dict[f"controls_{i}"] = r["controls"]
        sofa = r["sofa"]
        save_dict[f"boundary_{i}"] = (
            np.array(sofa.exterior.coords) if not sofa.is_empty else np.zeros((0, 2)))
        save_dict[f"corners_dense_{i}"] = r["corners_dense"]
        save_dict[f"eval_thetas_{i}"] = r["eval_thetas"]
    save_dict["areas"] = np.array(areas)
    np.savez(out_path, **save_dict)
    print(f"saved {out_path}")


if __name__ == "__main__":
    run()
