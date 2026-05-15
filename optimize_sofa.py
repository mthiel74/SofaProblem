"""Find a near-optimal moving-sofa shape.

Optimises the inner-corner trajectory C(theta) (M control points + cubic
spline) so that the intersection of corridor poses,
    S(C) = intersection over theta of L(theta; C(theta)),
has maximum area.
"""

from __future__ import annotations

import argparse
import time
import numpy as np
from scipy.optimize import minimize

from sofa import (
    build_sofa,
    control_thetas,
    evaluate_trajectory,
    hammersley_controls,
    neg_area_from_controls,
    quarter_arc_controls,
)


HAMMERSLEY_AREA = np.pi / 2.0 + 2.0 / np.pi      # ~ 2.2074
GERVER_AREA = 2.21953166887197                   # best known lower bound


def run(n_thetas: int, M: int, maxiter: int, seed: int, out_path: str,
        n_thetas_eval: int = 1600) -> None:
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, np.pi / 2.0, n_thetas)

    # Best seed in our experiments: corner sweeps a quarter circle of
    # radius 2/pi from (0, r) at theta = 0 to (r, 0) at theta = pi/2.
    init_controls = quarter_arc_controls(M, r=2.0 / np.pi)
    x0 = init_controls.reshape(-1).copy()

    best_x = x0
    best_neg = neg_area_from_controls(x0, thetas, M)
    print(f"M = {M} control points, n_thetas = {n_thetas} (eval at {n_thetas_eval}), vars = {x0.size}")
    print(f"initial area (optim) = {-best_neg:.6f}")

    # Two-stage: start from the seed and successively warm-start with the
    # best solution plus a small perturbation, to escape early local optima.
    t0 = time.time()
    current = x0
    for k in range(6):
        if k > 0:
            current = best_x + 0.05 * rng.standard_normal(best_x.size)
        res = minimize(
            neg_area_from_controls,
            current,
            args=(thetas, M),
            method="Nelder-Mead",
            options={
                "maxiter": maxiter,
                "xatol": 1e-6,
                "fatol": 1e-8,
                "adaptive": True,
                "disp": False,
            },
        )
        area = -res.fun
        elapsed = time.time() - t0
        improved = "*" if res.fun < best_neg else " "
        print(f" {improved} round {k}: area = {area:.6f}   (t = {elapsed:.1f}s, nfev = {res.nfev})")
        if res.fun < best_neg:
            best_neg = res.fun
            best_x = res.x

    best_area_optim = -best_neg
    # Re-evaluate the optimised trajectory at much denser theta sampling.
    # The optimisation samples are coarse, so the area at the optimisation
    # grid slightly overestimates the true sofa area (corridors that fall
    # between samples can still slice the sofa).
    control_xy = best_x.reshape(M, 2)
    thetas_eval = np.linspace(0.0, np.pi / 2.0, n_thetas_eval)
    corners_eval = evaluate_trajectory(control_xy, thetas_eval)
    sofa = build_sofa(corners_eval, thetas_eval)
    best_area_true = float(sofa.area)
    print()
    print(f"area at optimisation grid  : {best_area_optim:.6f}")
    print(f"area at dense eval grid    : {best_area_true:.6f}")
    print(f"Hammersley area            : {HAMMERSLEY_AREA:.6f}   (delta: {best_area_true - HAMMERSLEY_AREA:+.4f})")
    print(f"Gerver area (best known)   : {GERVER_AREA:.6f}   (delta: {best_area_true - GERVER_AREA:+.4f})")

    dense_thetas = np.linspace(0.0, np.pi / 2.0, 720)
    dense_corners = evaluate_trajectory(control_xy, dense_thetas)
    boundary = np.array(sofa.exterior.coords) if not sofa.is_empty else np.zeros((0, 2))
    holes = [np.array(r.coords) for r in (sofa.interiors if not sofa.is_empty else [])]

    np.savez(
        out_path,
        thetas=thetas,
        corners=evaluate_trajectory(control_xy, thetas),
        dense_thetas=dense_thetas,
        dense_corners=dense_corners,
        control_xy=control_xy,
        control_thetas=control_thetas(M),
        M=M,
        area=best_area_true,
        area_optim_grid=best_area_optim,
        boundary=boundary,
        holes=np.array(holes, dtype=object),
    )
    print(f"saved {out_path}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--thetas", type=int, default=100,
                   help="number of theta samples used while optimising")
    p.add_argument("--thetas-eval", type=int, default=1600,
                   help="number of theta samples used for the honest final-area readout")
    p.add_argument("--M", type=int, default=10,
                   help="number of cubic-spline control points on the corner trajectory")
    p.add_argument("--maxiter", type=int, default=3000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--out", type=str, default="sofa_result.npz")
    args = p.parse_args()
    run(args.thetas, args.M, args.maxiter, args.seed, args.out, args.thetas_eval)


if __name__ == "__main__":
    main()
