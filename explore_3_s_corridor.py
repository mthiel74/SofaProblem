"""Exploration 3: S-corridor (two consecutive 90 degree bends).

World corridor: arm 1 along +x at y in [0, 1], middle arm along +y at
x in [0, 1], arm 3 along -x at y in [M-1, M]. Two left turns, total
rotation pi. The sofa enters arm 1 going +x and exits arm 3 going -x.

We reuse the L-corridor framework with a wider 8-vertex polygon and a
theta range [0, pi] instead of [0, pi/2]. The trajectory of the FIRST
inner corner (1, 1) is what we optimise; the second inner corner sits
at (1, M) in world frame and is carried along by the rigid corridor
geometry.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from shapely.geometry import Polygon
from shapely.validation import make_valid


_BIG = 30.0
_M = 6.0  # middle-arm length; large enough not to bind the sofa

# World S-corridor, CCW
_WORLD_S = np.array([
    [0.0, 0.0],
    [_BIG, 0.0],
    [_BIG, 1.0],
    [1.0, 1.0],          # inner corner 1
    [1.0, _M],           # inner corner 2
    [-_BIG, _M],
    [-_BIG, _M - 1.0],
    [0.0, _M - 1.0],
], dtype=float)


def corridor_s(theta: float, cx: float, cy: float) -> Polygon:
    """S-corridor at orientation theta with inner corner 1 at (cx, cy)
    in the sofa-fixed frame. Same convention as the L-corridor code:
    world inner corner 1 (which is at (1, 1)) lands on (cx, cy)."""
    c, s = math.cos(theta), math.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    corner = np.array([cx, cy])
    pts = corner + (_WORLD_S[:, 0:1] - 1.0) * ex \
                 + (_WORLD_S[:, 1:2] - 1.0) * ey
    return Polygon(pts.tolist())


def build_sofa_s(corners: np.ndarray, thetas: np.ndarray) -> Polygon:
    sofa = corridor_s(float(thetas[0]), float(corners[0, 0]), float(corners[0, 1]))
    for i in range(1, len(thetas)):
        nxt = corridor_s(float(thetas[i]), float(corners[i, 0]), float(corners[i, 1]))
        sofa = sofa.intersection(nxt)
        if sofa.is_empty:
            return sofa
    return sofa if sofa.is_valid else make_valid(sofa)


def control_thetas(M: int, theta_max: float) -> np.ndarray:
    return np.linspace(0.0, theta_max, M)


def evaluate_trajectory(control_xy: np.ndarray, thetas: np.ndarray,
                         theta_max: float) -> np.ndarray:
    M = control_xy.shape[0]
    cth = control_thetas(M, theta_max)
    sx = CubicSpline(cth, control_xy[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cth, control_xy[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


def neg_area(flat: np.ndarray, thetas: np.ndarray, M: int, theta_max: float
             ) -> float:
    corners = evaluate_trajectory(flat.reshape(M, 2), thetas, theta_max)
    sofa = build_sofa_s(corners, thetas)
    return 0.0 if sofa.is_empty else -float(sofa.area)


def initial_trajectory(M: int, theta_max: float, r: float = 2.0 / math.pi
                        ) -> np.ndarray:
    """Quarter-arc-like seed extended to the longer theta range."""
    th = control_thetas(M, theta_max)
    # Wrap the quarter-arc twice — once per bend.
    phi = th * (math.pi / 2.0) / (theta_max / 2.0)  # scale so each "L stage" is one quarter
    return np.column_stack([r * np.sin(phi), r * np.cos(phi)])


def run(M: int = 14, n_thetas: int = 200, n_rounds: int = 4,
        maxiter: int = 3500, seed: int = 7,
        theta_max: float = math.pi) -> None:
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, theta_max, n_thetas)
    x0 = initial_trajectory(M, theta_max).reshape(-1)
    print(f"S-corridor: theta in [0, {math.degrees(theta_max):.0f} deg],  M = {M} ctrl pts")
    init_area = -neg_area(x0, thetas, M, theta_max)
    print(f"initial area = {init_area:.5f}")
    best, bf = x0, neg_area(x0, thetas, M, theta_max)
    for k in range(n_rounds):
        init = best if k == 0 else best + 0.05 * rng.standard_normal(best.size)
        res = minimize(neg_area, init, args=(thetas, M, theta_max),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-5, fatol=1e-7,
                                    adaptive=True))
        tag = "*" if res.fun < bf else " "
        print(f" {tag} round {k}: area = {-res.fun:.5f}")
        if res.fun < bf:
            bf, best = res.fun, res.x
    ctrl = best.reshape(M, 2)
    eval_thetas = np.linspace(0.0, theta_max, 1600)
    sofa = build_sofa_s(evaluate_trajectory(ctrl, eval_thetas, theta_max),
                        eval_thetas)
    print(f"\nfinal (honest)    area = {sofa.area:.5f}")
    print(f"Gerver single-bend area = 2.21953")
    print(f"S-corridor / single-bend ratio = {sofa.area / 2.21953:.4f}")
    boundary = np.array(sofa.exterior.coords) if not sofa.is_empty \
               else np.zeros((0, 2))
    np.savez("Wolfram/community-3d/data_s_corridor.npz",
             controls=ctrl, area=float(sofa.area),
             boundary=boundary, thetas=eval_thetas,
             corners=evaluate_trajectory(ctrl, eval_thetas, theta_max))
    print("saved data_s_corridor.npz")


if __name__ == "__main__":
    run()
