"""Optimise the moving sofa for a sequence of bend angles alpha.

Anchor angles chosen so they cover a wide span and give 7 numerical
results for the post table. Each optimisation takes 1-3 minutes.

Writes sweep_angles.npz containing, for each angle:
    alpha, area, controls (M x 2), sofa boundary, eval-grid trajectory.
"""

from __future__ import annotations

import math
import time
import numpy as np

from sofa_angle import optimise_sofa_alpha


ALPHAS_DEG = [60, 75, 90, 105, 120, 135, 150]


def run(out_path: str = "sweep_angles.npz") -> None:
    results = []
    t_start = time.time()
    warm = None
    for ad in ALPHAS_DEG:
        a = math.radians(ad)
        print(f"==== alpha = {ad} deg ====", flush=True)
        # Use the quarter-arc seed (warm-start across angles is unreliable
        # because theta_max changes with alpha and the controls'
        # interpretation shifts).
        res = optimise_sofa_alpha(
            alpha=a, M=10, n_thetas=80, n_rounds=3, maxiter=2200,
            warm_start=warm, verbose=True)
        elapsed = time.time() - t_start
        print(f"  -> area = {res['area']:.5f}    (cumulative t = {elapsed:.0f}s)", flush=True)
        results.append(res)

    # Bundle and save
    save_dict = {
        "alphas_deg": np.array([ad for ad in ALPHAS_DEG]),
        "areas": np.array([r["area"] for r in results]),
    }
    for i, r in enumerate(results):
        save_dict[f"controls_{i}"] = r["controls"]
        sofa = r["sofa"]
        if sofa.is_empty:
            save_dict[f"boundary_{i}"] = np.zeros((0, 2))
        else:
            save_dict[f"boundary_{i}"] = np.array(sofa.exterior.coords)
        save_dict[f"corners_dense_{i}"] = r["corners_dense"]
        save_dict[f"eval_thetas_{i}"] = r["eval_thetas"]
    np.savez(out_path, **save_dict)
    print(f"\nsaved {out_path}")


if __name__ == "__main__":
    run()
