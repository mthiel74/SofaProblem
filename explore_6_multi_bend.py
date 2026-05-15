"""Exploration 6: N consecutive 90 degree bends.

Build a corridor with N unit-width arms connected by N-1 left bends
(N = 1: classical L; N = 2: S-corridor; N = 3: a "C"; N = 4: U).
Optimise the inner-corner trajectory through theta in [0, N pi / 2].
Measure how the maximum sofa area S(N) scales with N.
"""

from __future__ import annotations

import math
import os
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from shapely.geometry import Polygon
from shapely.validation import make_valid


_BIG = 20.0


def world_polygon_nbend(N: int) -> np.ndarray:
    """N consecutive 90 deg LEFT bends; N - 1 inner corners on a corridor
    with N arms. Arm 1 in +x, then +y, then -x, then -y. Each arm
    extends _BIG units past its bend (so the arms don't bind the sofa)."""
    B = _BIG
    if N == 1:
        # L corridor, 1 bend, 2 arms
        return np.array([[0, 0], [B, 0], [B, 1], [1, 1], [1, B], [0, B]],
                        dtype=float)
    elif N == 2:
        # S corridor, 2 bends, 3 arms (+x, +y, -x)
        return np.array([
            [0, 0], [B, 0], [B, 1], [1, 1], [1, B],
            [-B, B], [-B, B - 1], [0, B - 1]
        ], dtype=float)
    elif N == 3:
        # C corridor, 3 bends, 4 arms (+x, +y, -x, -y)
        return np.array([
            [0, 0], [B, 0], [B, 1], [1, 1], [1, B],
            [-B, B], [-B, -B], [-B + 1, -B],
            [-B + 1, B - 1], [0, B - 1]
        ], dtype=float)
    elif N == 4:
        # Square loop, 4 bends, 5 arms (+x, +y, -x, -y, +x)
        return np.array([
            [0, 0], [B, 0], [B, 1], [1, 1], [1, B],
            [-B, B], [-B, -B], [-B + 1, -B],
            [-B + 1, -B + 1], [B, -B + 1],
            [B, -B + 2], [0, -B + 2]
            # tail closes the polygon back to (0, 0) via outer wall of arm 5
        ], dtype=float)
    raise NotImplementedError(f"N={N} not implemented")


def corridor_polygon(theta: float, cx: float, cy: float, N: int) -> Polygon:
    """N-bend corridor in sofa frame, world inner corner 1 at (cx, cy)."""
    c, s = math.cos(theta), math.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    WL = world_polygon_nbend(N)
    pts = np.array([cx, cy]) + (WL[:, 0:1] - 1) * ex + (WL[:, 1:2] - 1) * ey
    return Polygon(pts.tolist())


def build_sofa(corners: np.ndarray, thetas: np.ndarray, N: int) -> Polygon:
    sofa = corridor_polygon(float(thetas[0]), float(corners[0, 0]),
                             float(corners[0, 1]), N)
    for i in range(1, len(thetas)):
        nxt = corridor_polygon(float(thetas[i]), float(corners[i, 0]),
                                float(corners[i, 1]), N)
        sofa = sofa.intersection(nxt)
        if sofa.is_empty:
            return sofa
    return sofa if sofa.is_valid else make_valid(sofa)


def control_thetas(M: int, theta_max: float) -> np.ndarray:
    return np.linspace(0.0, theta_max, M)


def eval_traj(ctrl: np.ndarray, thetas: np.ndarray, theta_max: float) -> np.ndarray:
    M = ctrl.shape[0]
    cth = control_thetas(M, theta_max)
    sx = CubicSpline(cth, ctrl[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cth, ctrl[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


def neg_area(flat: np.ndarray, thetas: np.ndarray, M: int,
              N: int, theta_max: float) -> float:
    corners = eval_traj(flat.reshape(M, 2), thetas, theta_max)
    sofa = build_sofa(corners, thetas, N)
    return 0.0 if sofa.is_empty else -float(sofa.area)


def init_trajectory(M: int, N: int, theta_max: float,
                    r: float = 2.0 / math.pi) -> np.ndarray:
    """Wrap a quarter-arc N times across the longer theta range."""
    th = control_thetas(M, theta_max)
    phi = th * (math.pi / 2.0) / (theta_max / N)
    return np.column_stack([r * np.sin(phi), r * np.cos(phi)])


def run_one(N: int, M: int, n_thetas: int, n_rounds: int, maxiter: int,
            seed: int = 7) -> dict:
    theta_max = N * math.pi / 2.0
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, theta_max, n_thetas)
    x0 = init_trajectory(M, N, theta_max).reshape(-1)
    best, bf = x0, neg_area(x0, thetas, M, N, theta_max)
    print(f"  N={N}, M={M}, theta_max = {math.degrees(theta_max):.0f} deg")
    print(f"  initial area = {-bf:.5f}")
    for k in range(n_rounds):
        init = best if k == 0 else best + 0.05 * rng.standard_normal(best.size)
        res = minimize(neg_area, init, args=(thetas, M, N, theta_max),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-5, fatol=1e-7,
                                    adaptive=True))
        tag = "*" if res.fun < bf else " "
        print(f"   {tag} round {k}: area = {-res.fun:.5f}")
        if res.fun < bf:
            bf, best = res.fun, res.x
    ctrl = best.reshape(M, 2)
    eval_thetas = np.linspace(0.0, theta_max, 1200)
    sofa = build_sofa(eval_traj(ctrl, eval_thetas, theta_max), eval_thetas, N)
    return dict(N=N, area=float(sofa.area), controls=ctrl,
                sofa=sofa, theta_max=theta_max,
                eval_thetas=eval_thetas,
                corners_dense=eval_traj(ctrl, eval_thetas, theta_max))


def run():
    os.makedirs("Wolfram/community-3d", exist_ok=True)
    results = {}
    # N=1 baseline (use my existing Gerver-class number)
    results[1] = {"N": 1, "area": 2.21953,
                  "boundary": None,
                  "is_imported": True}
    print(f"==== N = 1 (reference) ====")
    print(f"  area = {results[1]['area']:.5f}  (Gerver, from prior exploration)")
    for N in [2, 3, 4]:
        print(f"==== N = {N} ====")
        # parameter scaling: more bends -> more control points
        M = 8 + 4 * (N - 1)        # 8, 12, 16
        n_thetas = 70 + 30 * (N - 1)  # 70, 100, 130
        n_rounds = 4
        maxiter = 2200
        r = run_one(N, M, n_thetas, n_rounds, maxiter)
        results[N] = r
        print(f"  -> S({N}) = {r['area']:.5f}")
    Ns = sorted(results.keys())
    Sns = [results[N]["area"] for N in Ns]
    print("\nN-bend scaling:")
    print(" N |   S(N)    | S(N)/S(1) | log(S(N)/S(1))")
    print("---|-----------|-----------|----------------")
    for N, A in zip(Ns, Sns):
        rel = A / Sns[0]
        print(f" {N} |  {A:.5f}  |   {rel:.4f}  |   {math.log(max(rel,1e-9)):.4f}")
    # save
    save = {"Ns": np.array(Ns), "areas": np.array(Sns)}
    for N in Ns:
        if not results[N].get("is_imported"):
            r = results[N]
            save[f"controls_{N}"] = r["controls"]
            sofa = r["sofa"]
            save[f"boundary_{N}"] = (np.array(sofa.exterior.coords)
                                     if not sofa.is_empty else np.zeros((0, 2)))
            save[f"corners_{N}"] = r["corners_dense"]
    np.savez("Wolfram/community-3d/data_multibend.npz", **save)
    print("\nsaved data_multibend.npz")


if __name__ == "__main__":
    run()
