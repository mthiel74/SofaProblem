"""Visualize the optimized sofa: static plot + animation of the corridor sweep."""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon as MplPolygon

from sofa import build_sofa, corridor_polygon, evaluate_trajectory


def _draw_sofa(ax, boundary: np.ndarray, holes, area: float) -> None:
    if boundary.size == 0:
        ax.set_title("(empty sofa)")
        return
    ax.add_patch(MplPolygon(boundary, closed=True, facecolor="#cfe2ff",
                            edgecolor="#1f4fa8", linewidth=2.0, zorder=2))
    for h in holes:
        ax.add_patch(MplPolygon(h, closed=True, facecolor="white",
                                edgecolor="#1f4fa8", linewidth=1.0, zorder=3))
    ax.set_aspect("equal")
    ax.set_title(f"optimized sofa, area = {area:.5f}")


def static_plot(data: dict, out_path: str) -> None:
    boundary = data["boundary"]
    holes = list(data["holes"])
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    _draw_sofa(ax, boundary, holes, float(data["area"]))

    corners = data["dense_corners"]
    ax.plot(corners[:, 0], corners[:, 1], color="#d2691e", linewidth=1.2,
            linestyle="--", label="inner-corner trajectory", zorder=4)
    ax.scatter([corners[0, 0], corners[-1, 0]],
               [corners[0, 1], corners[-1, 1]],
               color="#d2691e", s=25, zorder=5)
    ax.legend(loc="upper right", framealpha=0.9)

    # Comfortable padding around the sofa.
    xs = boundary[:, 0]
    ys = boundary[:, 1]
    pad = 0.5
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(ys.min() - pad, ys.max() + pad)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    print(f"saved {out_path}")
    plt.close(fig)


def animate(data: dict, out_path: str, n_frames: int = 90) -> None:
    """Animate the corridor sweeping over the sofa.

    For each theta, we draw the L-corridor (in the sofa-fixed frame) on top
    of the static sofa. The sofa is fully contained in every frame.
    """
    boundary = data["boundary"]
    holes = list(data["holes"])
    control_xy = data["control_xy"]

    frame_thetas = np.linspace(0.01, np.pi / 2 - 0.01, n_frames)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    xs = boundary[:, 0]
    ys = boundary[:, 1]
    pad = 1.2
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(ys.min() - pad, ys.max() + pad)
    ax.set_aspect("equal")
    sofa_patch = MplPolygon(boundary, closed=True, facecolor="#cfe2ff",
                            edgecolor="#1f4fa8", linewidth=2.0, zorder=2)
    ax.add_patch(sofa_patch)
    hole_patches = []
    for h in holes:
        hp = MplPolygon(h, closed=True, facecolor="white",
                        edgecolor="#1f4fa8", linewidth=1.0, zorder=3)
        ax.add_patch(hp)
        hole_patches.append(hp)

    corridor_patch = MplPolygon(np.zeros((4, 2)), closed=True,
                                facecolor="#fff0d6", edgecolor="#a04020",
                                linewidth=1.0, alpha=0.45, zorder=1)
    ax.add_patch(corridor_patch)
    title = ax.set_title("")

    def update(i):
        theta = frame_thetas[i]
        c = evaluate_trajectory(control_xy, np.array([theta]))[0]
        corr = corridor_polygon(theta, c[0], c[1])
        # Intersect with the visible window so the patch doesn't blow up.
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        from shapely.geometry import box as shp_box
        window = shp_box(xlim[0], ylim[0], xlim[1], ylim[1])
        clipped = corr.intersection(window)
        if clipped.is_empty:
            corridor_patch.set_xy(np.zeros((4, 2)))
        else:
            geom = clipped if clipped.geom_type == "Polygon" else max(
                clipped.geoms, key=lambda g: g.area)
            corridor_patch.set_xy(np.array(geom.exterior.coords))
        title.set_text(f"theta = {np.degrees(theta):5.1f} deg   "
                       f"area = {float(data['area']):.5f}")
        return corridor_patch, title

    anim = FuncAnimation(fig, update, frames=len(frame_thetas), interval=60, blit=False)
    writer = PillowWriter(fps=20)
    anim.save(out_path, writer=writer)
    print(f"saved {out_path}")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="sofa_result.npz")
    p.add_argument("--static-out", default="sofa.png")
    p.add_argument("--anim-out", default="sofa.gif")
    p.add_argument("--no-anim", action="store_true",
                   help="skip the animation (which is slower)")
    args = p.parse_args()

    npz = np.load(args.input, allow_pickle=True)
    data = {k: npz[k] for k in npz.files}
    static_plot(data, args.static_out)
    if not args.no_anim:
        animate(data, args.anim_out)


if __name__ == "__main__":
    main()
