"""Moving sofa with arbitrary corridor bend angle alpha in (0, pi).

Convention (interior bend angle, measured between the inner walls of the
two arms going outward from the inner corner):
  alpha = pi/2  ->  classical right-angle L
  alpha = pi    ->  straight strip (no bend)
  alpha -> 0    ->  folded back on itself

We keep the inner corner of the L fixed at (1, 1) in the world frame
(consistent with the right-angle code in sofa.py), so passing alpha=pi/2
must reproduce that polygon exactly. The trajectory of the inner corner
in the sofa-fixed frame, C(theta) = (X(theta), Y(theta)), is the
parameter we optimise.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from shapely.geometry import Polygon
from shapely.affinity import rotate as shp_rotate, translate as shp_translate
from shapely.validation import make_valid

_BIG = 40.0


def world_corridor_polygon(alpha: float) -> Polygon:
    """World-frame L-corridor with bend angle alpha, inner corner (1, 1),
    arm 1 along +x with floor y = 0 and ceiling y = 1.

    Geometry: inner walls leave the corner in directions u1 = (1, 0) and
    u2 = (cos alpha, sin alpha). Outer walls are parallel, 1 unit further
    away from the sofa side. The two outer walls meet at the outer
    corner O = (1 - (1 + cos alpha) / sin alpha, 0).
    """
    # Outer corner
    Ox = 1.0 - (1.0 + math.cos(alpha)) / math.sin(alpha)
    Oy = 0.0
    # Inner-arm-2 direction
    ca, sa = math.cos(alpha), math.sin(alpha)
    # Far ends along the two arms (parametrically _BIG units)
    inner_far_1 = (_BIG, 1.0)
    outer_far_1 = (_BIG, 0.0)
    inner_far_2 = (1.0 + _BIG * ca, 1.0 + _BIG * sa)
    outer_far_2 = (Ox + _BIG * ca, Oy + _BIG * sa)
    verts = [
        (Ox, Oy),       # outer corner
        outer_far_1,    # far end arm-1 floor
        inner_far_1,    # far end arm-1 ceiling
        (1.0, 1.0),     # inner corner
        inner_far_2,    # far end arm-2 inner wall
        outer_far_2,    # far end arm-2 outer wall
    ]
    poly = Polygon(verts)
    if not poly.is_valid:
        poly = make_valid(poly)
    return poly


def corridor_polygon_alpha(theta: float, cx: float, cy: float, alpha: float) -> Polygon:
    """L-corridor in the sofa-fixed frame, inner corner at (cx, cy)."""
    L = world_corridor_polygon(alpha)
    # In the sofa-fixed frame the world rotates by -theta around the world
    # corner (1, 1), then translates so (1, 1) lands on (cx, cy).
    Lr = shp_rotate(L, -math.degrees(theta), origin=(1.0, 1.0))
    return shp_translate(Lr, xoff=cx - 1.0, yoff=cy - 1.0)


def build_sofa_alpha(corners: np.ndarray, thetas: np.ndarray, alpha: float) -> Polygon:
    sofa = corridor_polygon_alpha(float(thetas[0]), float(corners[0, 0]), float(corners[0, 1]), alpha)
    for i in range(1, len(thetas)):
        nxt = corridor_polygon_alpha(float(thetas[i]), float(corners[i, 0]), float(corners[i, 1]), alpha)
        sofa = sofa.intersection(nxt)
        if sofa.is_empty:
            return sofa
    return sofa if sofa.is_valid else make_valid(sofa)


# ----------- trajectory parameterisation ----------------------------

def _control_thetas(M: int, theta_max: float) -> np.ndarray:
    return np.linspace(0.0, theta_max, M)


def evaluate_trajectory_alpha(control_xy: np.ndarray,
                              thetas: np.ndarray,
                              theta_max: float) -> np.ndarray:
    M = control_xy.shape[0]
    cthetas = _control_thetas(M, theta_max)
    sx = CubicSpline(cthetas, control_xy[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cthetas, control_xy[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


def neg_area_alpha(flat: np.ndarray, thetas: np.ndarray, M: int,
                    alpha: float, theta_max: float) -> float:
    corners = evaluate_trajectory_alpha(flat.reshape(M, 2), thetas, theta_max)
    sofa = build_sofa_alpha(corners, thetas, alpha)
    return 0.0 if sofa.is_empty else -float(sofa.area)


def quarter_arc_controls_alpha(M: int, theta_max: float, r: float = 2.0 / np.pi) -> np.ndarray:
    """Initial trajectory: corner sweeps an arc of radius r and angular
    span (pi/2) as theta runs from 0 to theta_max — this generalises the
    Hammersley seed of the right-angle case."""
    th = _control_thetas(M, theta_max)
    # Map theta in [0, theta_max] -> arc parameter in [0, pi/2] (so the
    # endpoints of the arc are always (0, r) and (r, 0), regardless of
    # theta_max).
    phi = th * (np.pi / 2.0) / theta_max if theta_max > 0 else th * 0.0
    return np.column_stack([r * np.sin(phi), r * np.cos(phi)])


def optimise_sofa_alpha(alpha: float,
                         M: int = 10,
                         n_thetas: int = 100,
                         n_rounds: int = 3,
                         maxiter: int = 2500,
                         seed: int = 7,
                         warm_start: np.ndarray | None = None,
                         verbose: bool = True) -> dict:
    """Optimise the corner trajectory for bend angle alpha.

    Returns a dict with controls (Mx2), area (high-resolution), the
    sofa polygon, and dense-sampled corner trajectory.
    """
    rng = np.random.default_rng(seed)
    theta_max = np.pi - alpha
    thetas = np.linspace(0.0, theta_max, n_thetas)
    x0 = (warm_start if warm_start is not None
          else quarter_arc_controls_alpha(M, theta_max)).reshape(-1).copy()

    best_x, best_neg = x0, neg_area_alpha(x0, thetas, M, alpha, theta_max)
    if verbose:
        print(f"  alpha = {math.degrees(alpha):5.1f} deg, theta_max = {math.degrees(theta_max):5.1f} deg")
        print(f"  initial area = {-best_neg:.5f}")
    for k in range(n_rounds):
        init = best_x if k == 0 else best_x + 0.05 * rng.standard_normal(best_x.size)
        res = minimize(neg_area_alpha, init,
                       args=(thetas, M, alpha, theta_max),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-5, fatol=1e-7, adaptive=True))
        tag = "*" if res.fun < best_neg else " "
        if verbose:
            print(f"  {tag} round {k}: area = {-res.fun:.5f}")
        if res.fun < best_neg:
            best_neg, best_x = res.fun, res.x

    control_xy = best_x.reshape(M, 2)
    eval_thetas = np.linspace(0.0, theta_max, 1200)
    sofa = build_sofa_alpha(
        evaluate_trajectory_alpha(control_xy, eval_thetas, theta_max),
        eval_thetas, alpha)
    return dict(alpha=alpha, theta_max=theta_max,
                controls=control_xy, area=float(sofa.area),
                sofa=sofa, eval_thetas=eval_thetas,
                corners_dense=evaluate_trajectory_alpha(control_xy, eval_thetas, theta_max))
