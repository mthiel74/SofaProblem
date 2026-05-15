"""End-to-end moving-sofa pipeline (self-contained).

Run as a single script:
    python3 sofa_pipeline.py
to produce sofa.png and sofa.gif. Reproduces the result reported in the
Wolfram Community post: an optimised "Gerver-class" sofa of area ≈ 2.2197
(Gerver's lower bound is 2.21953).

This is the same code that lives in sofa.py + optimize_sofa.py +
visualize.py, but pulled into one file so it can be dropped into a
notebook cell or a README.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon as MplPolygon
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from shapely.geometry import Polygon, box as shp_box
from shapely.validation import make_valid

# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------

_BIG = 40.0  # bounding box for the (otherwise infinite) corridor arms

# World L-corridor: first quadrant minus the {x>1, y>1} block.
_WORLD_L = np.array([
    [0.0, 0.0], [_BIG, 0.0], [_BIG, 1.0],
    [1.0, 1.0], [1.0, _BIG], [0.0, _BIG],
])


def corridor_polygon(theta, cx, cy):
    """L-corridor at orientation theta, inner corner at (cx, cy) in sofa frame.

    In the sofa-fixed frame, world x-axis -> e_x(theta) = (cos t, -sin t),
    world y-axis -> e_y(theta) = (sin t, cos t). A world point (xw, yw)
    maps to corner + (xw-1) * e_x + (yw-1) * e_y, so the world inner
    corner (1, 1) lands on (cx, cy).
    """
    c, s = np.cos(theta), np.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    corner = np.array([cx, cy])
    pts = corner + (_WORLD_L[:, 0:1] - 1.0) * ex + (_WORLD_L[:, 1:2] - 1.0) * ey
    return Polygon(pts.tolist())


def build_sofa(corners, thetas):
    """Intersect every corridor pose -> the largest sofa for this trajectory."""
    sofa = corridor_polygon(float(thetas[0]), float(corners[0, 0]), float(corners[0, 1]))
    for i in range(1, len(thetas)):
        sofa = sofa.intersection(
            corridor_polygon(float(thetas[i]), float(corners[i, 0]), float(corners[i, 1])))
        if sofa.is_empty:
            return sofa
    return sofa if sofa.is_valid else make_valid(sofa)


# ----------------------------------------------------------------------
# Trajectory: M control points -> cubic spline
# ----------------------------------------------------------------------

def control_thetas(M):
    return np.linspace(0.0, np.pi / 2.0, M)

def evaluate_trajectory(control_xy, thetas):
    M = control_xy.shape[0]
    cthetas = control_thetas(M)
    sx = CubicSpline(cthetas, control_xy[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cthetas, control_xy[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])

def neg_area_from_controls(flat, thetas, M):
    corners = evaluate_trajectory(flat.reshape(M, 2), thetas)
    sofa = build_sofa(corners, thetas)
    return 0.0 if sofa.is_empty else -float(sofa.area)


def quarter_arc_controls(M, r=2.0 / np.pi):
    """Initial trajectory: corner sweeps a quarter circle of radius r."""
    th = control_thetas(M)
    return np.column_stack([r * np.sin(th), r * np.cos(th)])


# ----------------------------------------------------------------------
# Optimisation
# ----------------------------------------------------------------------

def optimise_sofa(M=10, n_thetas=100, n_rounds=4, maxiter=3000, seed=7):
    rng = np.random.default_rng(seed)
    thetas = np.linspace(0.0, np.pi / 2.0, n_thetas)
    x0 = quarter_arc_controls(M).reshape(-1)
    best_x, best_neg = x0, neg_area_from_controls(x0, thetas, M)
    print(f"initial area = {-best_neg:.5f}")
    for k in range(n_rounds):
        init = best_x if k == 0 else best_x + 0.05 * rng.standard_normal(best_x.size)
        res = minimize(neg_area_from_controls, init, args=(thetas, M),
                       method="Nelder-Mead",
                       options=dict(maxiter=maxiter, xatol=1e-5, fatol=1e-7, adaptive=True))
        tag = "*" if res.fun < best_neg else " "
        print(f" {tag} round {k}: area = {-res.fun:.5f}")
        if res.fun < best_neg:
            best_neg, best_x = res.fun, res.x
    return best_x.reshape(M, 2)


# ----------------------------------------------------------------------
# Visualisation
# ----------------------------------------------------------------------

def static_plot(control_xy, out_path="sofa.png"):
    dense_thetas = np.linspace(0.0, np.pi / 2.0, 720)
    eval_thetas = np.linspace(0.0, np.pi / 2.0, 1600)  # honest area
    corners = evaluate_trajectory(control_xy, dense_thetas)
    sofa = build_sofa(evaluate_trajectory(control_xy, eval_thetas), eval_thetas)
    boundary = np.array(sofa.exterior.coords)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.add_patch(MplPolygon(boundary, closed=True, facecolor="#cfe2ff",
                            edgecolor="#1f4fa8", linewidth=2.0, zorder=2))
    for r in sofa.interiors:
        ax.add_patch(MplPolygon(np.array(r.coords), closed=True, facecolor="white",
                                edgecolor="#1f4fa8", linewidth=1.0, zorder=3))
    ax.plot(corners[:, 0], corners[:, 1], color="#d2691e",
            linewidth=1.2, linestyle="--", label="inner-corner trajectory")
    ax.scatter([corners[0, 0], corners[-1, 0]], [corners[0, 1], corners[-1, 1]],
               color="#d2691e", s=25)
    ax.set_aspect("equal"); ax.legend(loc="upper right", framealpha=0.9)
    xs, ys = boundary[:, 0], boundary[:, 1]
    ax.set_xlim(xs.min() - 0.5, xs.max() + 0.5)
    ax.set_ylim(ys.min() - 0.5, ys.max() + 0.5)
    ax.set_title(f"optimised sofa, area = {sofa.area:.5f}")
    fig.tight_layout(); fig.savefig(out_path, dpi=160); plt.close(fig)
    return float(sofa.area), boundary, list(sofa.interiors)


def animate(control_xy, out_path="sofa.gif", n_frames=90):
    eval_thetas = np.linspace(0.0, np.pi / 2.0, 1600)
    sofa = build_sofa(evaluate_trajectory(control_xy, eval_thetas), eval_thetas)
    boundary = np.array(sofa.exterior.coords)
    holes = [np.array(r.coords) for r in sofa.interiors]
    frame_thetas = np.linspace(0.01, np.pi / 2 - 0.01, n_frames)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    xs, ys = boundary[:, 0], boundary[:, 1]
    ax.set_xlim(xs.min() - 1.2, xs.max() + 1.2)
    ax.set_ylim(ys.min() - 1.2, ys.max() + 1.2)
    ax.set_aspect("equal")
    ax.add_patch(MplPolygon(boundary, closed=True, facecolor="#cfe2ff",
                            edgecolor="#1f4fa8", linewidth=2.0, zorder=2))
    for h in holes:
        ax.add_patch(MplPolygon(h, closed=True, facecolor="white",
                                edgecolor="#1f4fa8", linewidth=1.0, zorder=3))
    corridor_patch = MplPolygon(np.zeros((4, 2)), closed=True,
                                facecolor="#fff0d6", edgecolor="#a04020",
                                linewidth=1.0, alpha=0.45, zorder=1)
    ax.add_patch(corridor_patch)
    title = ax.set_title("")
    window = shp_box(*ax.get_xlim(), *ax.get_ylim())  # type: ignore

    def update(i):
        theta = frame_thetas[i]
        c = evaluate_trajectory(control_xy, np.array([theta]))[0]
        clipped = corridor_polygon(theta, c[0], c[1]).intersection(window)
        if clipped.is_empty:
            corridor_patch.set_xy(np.zeros((4, 2)))
        else:
            geom = clipped if clipped.geom_type == "Polygon" else \
                max(clipped.geoms, key=lambda g: g.area)
            corridor_patch.set_xy(np.array(geom.exterior.coords))
        title.set_text(f"theta = {np.degrees(theta):5.1f}°   area = {sofa.area:.5f}")
        return corridor_patch, title

    anim = FuncAnimation(fig, update, frames=len(frame_thetas), interval=60, blit=False)
    anim.save(out_path, writer=PillowWriter(fps=20))
    plt.close(fig)


# ----------------------------------------------------------------------

if __name__ == "__main__":
    control_xy = optimise_sofa(M=10, n_thetas=100, n_rounds=4)
    area, _, _ = static_plot(control_xy, "sofa.png")
    animate(control_xy, "sofa.gif")
    print(f"\nfinal area (1600-sample eval): {area:.5f}")
    print(f"Hammersley (analytic)        : {np.pi/2 + 2/np.pi:.5f}")
    print(f"Gerver (best known)          : 2.21953")
