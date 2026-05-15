"""Exploration 5: Largest CONVEX moving sofa.

For convex sofas only, the half-disk of radius 1 is known to be
optimal (it inscribes in every corridor pose and area = pi/2).
We reproduce this numerically by adding a convexity-violating
penalty to the optimiser. Convergence should land at pi/2.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize
from shapely.geometry import Polygon

from sofa import build_sofa, evaluate_trajectory, control_thetas


def convexity_violation(boundary: np.ndarray) -> float:
    """Sum of negative cross products around the polygon. Zero iff convex."""
    n = len(boundary) - 1  # last point repeats the first
    violation = 0.0
    for i in range(n):
        a = boundary[i]
        b = boundary[(i + 1) % n]
        c = boundary[(i + 2) % n]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if cross < 0:
            violation += -cross
    return violation


def neg_area_convex(flat: np.ndarray, thetas: np.ndarray, M: int,
                    lam: float = 5.0) -> float:
    """Negative sofa area with a penalty for any non-convex boundary."""
    corners = evaluate_trajectory(flat.reshape(M, 2), thetas)
    sofa = build_sofa(corners, thetas)
    if sofa.is_empty:
        return 0.0
    if sofa.geom_type == "MultiPolygon":
        # take the largest connected piece
        sofa = max(sofa.geoms, key=lambda g: g.area)
    boundary = np.array(sofa.exterior.coords)
    penalty = convexity_violation(boundary)
    return -float(sofa.area) + lam * penalty


def diameter_sweep_init(M: int) -> np.ndarray:
    """The inner corner of the half-disk sofa traces the disk's diameter
    as the sofa moves through the L. Initialise the trajectory to be
    that straight line from (1, 0) at theta = 0 to (-1, 0) at theta = pi/2.
    """
    th = control_thetas(M)
    return np.column_stack([1.0 - 2.0 * th / np.pi, np.zeros_like(th)])


def quarter_arc_init(M: int, r: float = 0.5) -> np.ndarray:
    """Small-radius quarter arc fallback."""
    th = control_thetas(M)
    return np.column_stack([r * np.sin(th), r * np.cos(th)])


def run(M: int = 8, n_thetas: int = 100, n_rounds: int = 4,
        maxiter: int = 3000, seed: int = 13):
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, np.pi / 2.0, n_thetas)
    x0 = diameter_sweep_init(M).reshape(-1)
    print(f"target (Gibbs 1968 half-disk) : pi/2 = {math.pi / 2:.6f}")
    print(f"initial : area_penalised = {-neg_area_convex(x0, thetas, M):.5f}")
    best, bf = x0, neg_area_convex(x0, thetas, M)
    for k in range(n_rounds):
        init = best if k == 0 else best + 0.04 * rng.standard_normal(best.size)
        res = minimize(neg_area_convex, init, args=(thetas, M),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-6, fatol=1e-8,
                                    adaptive=True))
        tag = "*" if res.fun < bf else " "
        print(f" {tag} round {k}: penalised = {-res.fun:.5f}")
        if res.fun < bf:
            bf, best = res.fun, res.x
    ctrl = best.reshape(M, 2)
    # Honest area (no penalty) of the resulting sofa
    sofa = build_sofa(evaluate_trajectory(ctrl, np.linspace(0, np.pi / 2, 1600)),
                      np.linspace(0, np.pi / 2, 1600))
    if sofa.geom_type == "MultiPolygon":
        sofa = max(sofa.geoms, key=lambda g: g.area)
    print(f"\nfinal area (honest) = {sofa.area:.5f}")
    print(f"target              = {math.pi / 2:.5f}")
    boundary = np.array(sofa.exterior.coords)
    print(f"convex violation    = {convexity_violation(boundary):.6f}")
    np.savez("Wolfram/community-3d/data_convex.npz",
             controls=ctrl, area=float(sofa.area),
             boundary=boundary)
    print("saved data_convex.npz")


if __name__ == "__main__":
    run()
