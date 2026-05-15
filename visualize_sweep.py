"""Build the angle-sweep figures from sweep_angles.npz:
    sweep_static.png  - all 7 anchor sofas in a row
    sweep.gif         - animation across alpha (linear interp of controls)
"""

from __future__ import annotations

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon as MplPolygon

from sofa_angle import (
    build_sofa_alpha,
    evaluate_trajectory_alpha,
    world_corridor_polygon,
)


def _draw_sofa_and_corridor(ax, sofa_boundary, alpha, title):
    """Draw the corridor (clipped to a window around the sofa) and the sofa."""
    from shapely.geometry import box as shp_box
    ax.cla()
    ax.set_aspect("equal")

    if sofa_boundary.size:
        xs, ys = sofa_boundary[:, 0], sofa_boundary[:, 1]
        xmin, xmax, ymin, ymax = xs.min(), xs.max(), ys.min(), ys.max()
    else:
        xmin, xmax, ymin, ymax = -1.0, 2.0, -1.0, 2.0
    pad = max(0.6, 0.15 * max(xmax - xmin, ymax - ymin))
    xlo, xhi = xmin - pad, xmax + pad
    ylo, yhi = ymin - pad, ymax + pad

    L = world_corridor_polygon(alpha)
    L_clip = L.intersection(shp_box(xlo, ylo, xhi, yhi))
    if not L_clip.is_empty:
        geoms = [L_clip] if L_clip.geom_type == "Polygon" else list(L_clip.geoms)
        for g in geoms:
            cx, cy = g.exterior.xy
            ax.fill(cx, cy, color="#f6f6f6", edgecolor="#999999", linewidth=0.8, zorder=0)

    if sofa_boundary.size:
        ax.add_patch(MplPolygon(sofa_boundary, closed=True, facecolor="#cfe2ff",
                                edgecolor="#1f4fa8", linewidth=1.6, zorder=2))
    ax.set_title(title, fontsize=10)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.tick_params(labelsize=7)


def static_sweep(npz_path: str = "sweep_angles.npz",
                 out_path: str = "sweep_static.png") -> None:
    data = np.load(npz_path, allow_pickle=True)
    angles = data["alphas_deg"]
    areas = data["areas"]
    n = len(angles)
    fig, axes = plt.subplots(1, n, figsize=(2.4 * n, 3.4))
    if n == 1:
        axes = [axes]
    for i, (ax, ad, a) in enumerate(zip(axes, angles, areas)):
        boundary = data[f"boundary_{i}"]
        _draw_sofa_and_corridor(ax, boundary, math.radians(float(ad)),
                                f"α = {int(ad)}°\narea = {a:.3f}")
    fig.suptitle("Optimised sofa at varying bend angle α", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


def sweep_animation(npz_path: str = "sweep_angles.npz",
                    out_path: str = "sweep.gif",
                    n_frames: int = 90) -> None:
    """Animate sofa as alpha varies. Anchors are the optimised sofas;
    between them the control points are linearly interpolated."""
    data = np.load(npz_path, allow_pickle=True)
    angles = data["alphas_deg"].astype(float)
    n_anchors = len(angles)
    controls = [data[f"controls_{i}"] for i in range(n_anchors)]
    M = controls[0].shape[0]

    # frame angles span the anchor range
    alphas_frames = np.linspace(angles[0], angles[-1], n_frames)

    def interpolate_controls(alpha_deg: float) -> np.ndarray:
        # piecewise linear in the anchor sequence
        if alpha_deg <= angles[0]:
            return controls[0]
        if alpha_deg >= angles[-1]:
            return controls[-1]
        i = np.searchsorted(angles, alpha_deg) - 1
        t = (alpha_deg - angles[i]) / (angles[i + 1] - angles[i])
        return (1 - t) * controls[i] + t * controls[i + 1]

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    title_left = fig.text(0.5, 0.95, "", ha="center", fontsize=11)

    def update(idx):
        ad = float(alphas_frames[idx])
        a = math.radians(ad)
        theta_max = math.pi - a
        ctrl = interpolate_controls(ad)
        eval_thetas = np.linspace(0.0, theta_max, 240)
        corners = evaluate_trajectory_alpha(ctrl, eval_thetas, theta_max)
        sofa = build_sofa_alpha(corners, eval_thetas, a)
        boundary = np.array(sofa.exterior.coords) if (not sofa.is_empty) else np.zeros((0, 2))
        area = float(sofa.area) if not sofa.is_empty else 0.0
        _draw_sofa_and_corridor(ax, boundary, a,
                                f"α = {ad:5.1f}°    area = {area:.3f}")
        return ax,

    anim = FuncAnimation(fig, update, frames=n_frames, interval=80, blit=False)
    anim.save(out_path, writer=PillowWriter(fps=15))
    plt.close(fig)
    print(f"saved {out_path}")


if __name__ == "__main__":
    static_sweep()
    sweep_animation()
