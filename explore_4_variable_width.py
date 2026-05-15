"""Exploration 4: Asymmetric-width L-corridor.

Arm 1 has width 1, arm 2 has width w in (0, 1]. Inner corner at (1, 1)
in world frame. World L:
  Arm 1: y in [0, 1],   x in [1, BIG]
  Arm 2: x in [1-w, 1], y in [1, BIG]
  Outer corner: (1-w, 0)

Recovers Gerver at w = 1, shrinks as w -> 0. Same framework as the
right-angle sofa optimiser; only the world polygon and the inner-corner
trajectory change.
"""

from __future__ import annotations

import math
import os
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from shapely.geometry import Polygon
from shapely.validation import make_valid


_BIG = 40.0


def world_l_w(w: float) -> np.ndarray:
    return np.array([
        [1.0 - w, 0.0], [_BIG, 0.0], [_BIG, 1.0],
        [1.0, 1.0], [1.0, _BIG], [1.0 - w, _BIG],
    ], dtype=float)


def corridor_w(theta: float, cx: float, cy: float, w: float) -> Polygon:
    c, s = math.cos(theta), math.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    corner = np.array([cx, cy])
    WL = world_l_w(w)
    pts = corner + (WL[:, 0:1] - 1.0) * ex + (WL[:, 1:2] - 1.0) * ey
    return Polygon(pts.tolist())


def build_sofa_w(corners: np.ndarray, thetas: np.ndarray, w: float) -> Polygon:
    sofa = corridor_w(float(thetas[0]), float(corners[0, 0]), float(corners[0, 1]), w)
    for i in range(1, len(thetas)):
        nxt = corridor_w(float(thetas[i]), float(corners[i, 0]), float(corners[i, 1]), w)
        sofa = sofa.intersection(nxt)
        if sofa.is_empty:
            return sofa
    return sofa if sofa.is_valid else make_valid(sofa)


def control_thetas(M: int) -> np.ndarray:
    return np.linspace(0.0, math.pi / 2.0, M)


def eval_traj(ctrl: np.ndarray, thetas: np.ndarray) -> np.ndarray:
    M = ctrl.shape[0]
    cth = control_thetas(M)
    sx = CubicSpline(cth, ctrl[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cth, ctrl[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


def neg_area(flat: np.ndarray, thetas: np.ndarray, M: int, w: float) -> float:
    corners = eval_traj(flat.reshape(M, 2), thetas)
    sofa = build_sofa_w(corners, thetas, w)
    return 0.0 if sofa.is_empty else -float(sofa.area)


def quarter_arc(M: int, r: float = 2.0 / math.pi) -> np.ndarray:
    th = control_thetas(M)
    return np.column_stack([r * np.sin(th), r * np.cos(th)])


def optimise(w: float, M: int = 10, n_thetas: int = 100, n_rounds: int = 4,
             maxiter: int = 2500, seed: int = 7,
             warm_start: np.ndarray | None = None,
             verbose: bool = True) -> dict:
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, math.pi / 2.0, n_thetas)
    x0 = warm_start.reshape(-1) if warm_start is not None else quarter_arc(M).reshape(-1)
    best, bf = x0, neg_area(x0, thetas, M, w)
    if verbose:
        print(f"  w = {w:.3f}, initial area = {-bf:.5f}")
    for k in range(n_rounds):
        init = best if k == 0 else best + 0.05 * rng.standard_normal(best.size)
        res = minimize(neg_area, init, args=(thetas, M, w),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-5, fatol=1e-7,
                                    adaptive=True))
        tag = "*" if res.fun < bf else " "
        if verbose:
            print(f"    {tag} round {k}: area = {-res.fun:.5f}")
        if res.fun < bf:
            bf, best = res.fun, res.x
    ctrl = best.reshape(M, 2)
    eval_thetas = np.linspace(0.0, math.pi / 2.0, 1200)
    sofa = build_sofa_w(eval_traj(ctrl, eval_thetas), eval_thetas, w)
    return dict(w=w, area=float(sofa.area), controls=ctrl,
                sofa=sofa, eval_thetas=eval_thetas,
                corners_dense=eval_traj(ctrl, eval_thetas))


def run():
    os.makedirs("Wolfram/community-3d", exist_ok=True)
    widths = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]  # cascade from w = 1 down
    results = {}
    warm = None
    for w in widths:
        print(f"==== width w = {w} ====")
        r = optimise(w, M=10, n_thetas=80, n_rounds=3, maxiter=2000,
                     warm_start=warm, verbose=True)
        warm = r["controls"]
        print(f"  -> area(w={w}) = {r['area']:.5f}")
        results[w] = r
    save = {"widths": np.array(widths),
            "areas": np.array([results[w]["area"] for w in widths])}
    for i, w in enumerate(widths):
        save[f"controls_{i}"] = results[w]["controls"]
        save[f"boundary_{i}"] = (
            np.array(results[w]["sofa"].exterior.coords)
            if not results[w]["sofa"].is_empty else np.zeros((0, 2)))
    np.savez("Wolfram/community-3d/data_variable_width.npz", **save)
    print("\nfinal table:")
    print("  w      area")
    for w in widths:
        print(f"  {w:.2f}   {results[w]['area']:.5f}")
    print("\nsaved data_variable_width.npz")


if __name__ == "__main__":
    run()
