"""Moving-sofa problem: find the largest 2D shape that can navigate an L-corridor.

Formulation (sofa-fixed frame):
  Instead of moving the sofa through a fixed L-shaped corridor, we fix the
  sofa and let the L-corridor rotate/translate around it. The motion is
  parameterized by an angle theta in [0, pi/2]. In the sofa-fixed frame,
  the inner corner of the L follows a trajectory
        C(theta) = (X(theta), Y(theta)).
  For a chosen trajectory, the largest sofa that survives the entire motion
  is the polygon intersection of every L-shaped corridor pose:
        S*(C) = intersection over theta of  L(theta; C(theta)).
  We then maximize  area(S*(C))  over the trajectory C.

  The world-frame L-corridor is the standard one:
        {(x, y) : x >= 0, y >= 0,  NOT (x > 1 AND y > 1)}
  i.e. the first quadrant with the {x>1, y>1} corner cut out. Width 1.

  The trajectory C is represented by a small number of control points
  in (X, Y), connected by a periodic-free cubic interpolation. The
  optimizer adjusts those control points.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline
from shapely.geometry import Polygon
from shapely.validation import make_valid

# Bounding box for the (otherwise infinite) corridor arms. Plenty larger than
# any realistic sofa, which lives within a few units of the corner.
_BIG = 40.0

# Vertices of the world-frame L-corridor (counter-clockwise), clipped to _BIG.
# Inner corner at (1, 1).
_WORLD_L_VERTICES = np.array([
    [0.0, 0.0],
    [_BIG, 0.0],
    [_BIG, 1.0],
    [1.0, 1.0],
    [1.0, _BIG],
    [0.0, _BIG],
])


def corridor_polygon(theta: float, cx: float, cy: float) -> Polygon:
    """L-corridor at parameter theta with inner corner at (cx, cy) in sofa frame.

    In the sofa-fixed frame, the world's x-axis points along
        e_x(theta) = (cos theta, -sin theta)
    and the world's y-axis points along
        e_y(theta) = (sin theta,  cos theta).
    A world-frame point (xw, yw) maps to
        corner + (xw - 1) * e_x + (yw - 1) * e_y
    so that world point (1, 1) (the inner corner) lands on (cx, cy).
    """
    c, s = np.cos(theta), np.sin(theta)
    ex = np.array([c, -s])
    ey = np.array([s, c])
    corner = np.array([cx, cy])
    pts = corner + (_WORLD_L_VERTICES[:, 0:1] - 1.0) * ex \
                 + (_WORLD_L_VERTICES[:, 1:2] - 1.0) * ey
    return Polygon(pts.tolist())


def build_sofa(corners: np.ndarray, thetas: np.ndarray) -> Polygon:
    """Intersect every corridor pose. corners has shape (N, 2), one per theta."""
    if len(thetas) == 0:
        return Polygon()
    sofa = corridor_polygon(float(thetas[0]), float(corners[0, 0]), float(corners[0, 1]))
    for i in range(1, len(thetas)):
        nxt = corridor_polygon(float(thetas[i]), float(corners[i, 0]), float(corners[i, 1]))
        sofa = sofa.intersection(nxt)
        if sofa.is_empty:
            return sofa
    if not sofa.is_valid:
        sofa = make_valid(sofa)
    return sofa


# ---------------------------------------------------------------------------
# Trajectory parameterization: M control points (X_k, Y_k) at fixed theta
# locations, connected by a cubic spline. We optimize the 2*M (X, Y) values.
#
# Note: shifting C(theta) by a constant vector v translates the sofa by v
# without changing its area, so there are 2 gauge degrees of freedom in the
# parameterization. The optimizer will wander along them without consequence;
# we simply post-translate the optimised polygon if a canonical placement is
# wanted.
# ---------------------------------------------------------------------------


def control_thetas(M: int) -> np.ndarray:
    """Theta values of the M control points, evenly spaced over [0, pi/2]."""
    return np.linspace(0.0, np.pi / 2.0, M)


def evaluate_trajectory(control_xy: np.ndarray, thetas: np.ndarray) -> np.ndarray:
    """Cubic-spline interpolate the control points to the requested thetas.

    control_xy has shape (M, 2): the (X, Y) values at the M control thetas.
    Returns (len(thetas), 2).
    """
    M = control_xy.shape[0]
    cthetas = control_thetas(M)
    sx = CubicSpline(cthetas, control_xy[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cthetas, control_xy[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


def neg_area_from_controls(flat: np.ndarray,
                           thetas: np.ndarray,
                           M: int) -> float:
    """Objective for scipy.optimize: -area(sofa). flat has 2M floats."""
    control_xy = flat.reshape(M, 2)
    corners = evaluate_trajectory(control_xy, thetas)
    sofa = build_sofa(corners, thetas)
    if sofa.is_empty:
        return 0.0
    return -float(sofa.area)


def quarter_arc_controls(M: int, r: float = 2.0 / np.pi) -> np.ndarray:
    """Simplest seed: corner sweeps a quarter circle of radius r,
        X(theta) = r * sin(theta),   Y(theta) = r * cos(theta),
    going from (0, r) at theta = 0 to (r, 0) at theta = pi/2.
    """
    thetas = control_thetas(M)
    X = r * np.sin(thetas)
    Y = r * np.cos(thetas)
    return np.column_stack([X, Y])


def hammersley_controls(M: int, x_far: float = -2.5) -> np.ndarray:
    """Hammersley-style seed with the corner pulled far to -x at the
    endpoints so the L's forbidden-corner block stays well clear of the
    sofa, with the corner approaching x = -r, y = 1/2 at theta = pi/4.
    """
    r = 2.0 / np.pi
    thetas = control_thetas(M)
    X = -r + (x_far + r) * (1.0 - np.sin(2.0 * thetas))
    Y = np.cos(thetas) ** 2
    return np.column_stack([X, Y])
