"""Build the per-angle aisle/sofa thumbnails for the post's Section-7
table, an area-vs-angle plot, and a 3D extruded-sofa visualization for
the higher-dimensions section.

Outputs (Wolfram/community/):
  aisle_<deg>.png        small corridor outline, one per angle
  sofa_<deg>.png         small sofa shape, one per angle
  area_vs_angle.png      curve of area S(alpha)
  gerver_3d.png          extruded 2D Gerver sofa as a 3D surface
  thumb_unit_square.png  thumbnails for Table 1
  thumb_half_disk.png
  thumb_hammersley.png
  thumb_gerver.png       (same as the optimised sofa)
"""

from __future__ import annotations

import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, FancyBboxPatch, Wedge
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from shapely.geometry import box as shp_box

from sofa_angle import world_corridor_polygon


OUTDIR = "Wolfram/community"
os.makedirs(OUTDIR, exist_ok=True)


def _save(fig, name):
    fig.savefig(os.path.join(OUTDIR, name), dpi=120, bbox_inches="tight",
                facecolor="white", transparent=False)
    plt.close(fig)


def _scale_for_all_sofas(sweep_data):
    """Common bounding box for all sofas + aisles, so the thumbnails are
    visually comparable."""
    n = len(sweep_data["alphas_deg"])
    xmin = ymin = float("inf")
    xmax = ymax = float("-inf")
    for i in range(n):
        b = sweep_data[f"boundary_{i}"]
        if b.size == 0:
            continue
        xmin = min(xmin, b[:, 0].min())
        xmax = max(xmax, b[:, 0].max())
        ymin = min(ymin, b[:, 1].min())
        ymax = max(ymax, b[:, 1].max())
    return xmin, xmax, ymin, ymax


def aisle_and_sofa_thumbs(npz_path: str = "sweep_angles.npz") -> None:
    data = np.load(npz_path, allow_pickle=True)
    angles = data["alphas_deg"].astype(int)
    xmin, xmax, ymin, ymax = _scale_for_all_sofas(data)
    pad = 0.6
    # use the same aspect across all images
    span_x = xmax - xmin + 2 * pad
    span_y = ymax - ymin + 2 * pad

    for i, ad in enumerate(angles):
        a = math.radians(float(ad))

        # ----- aisle ------
        figA, axA = plt.subplots(figsize=(2.6, 2.6))
        axA.set_aspect("equal")
        axA.set_xticks([]); axA.set_yticks([])
        for s in axA.spines.values():
            s.set_visible(False)
        # Use a smaller drawing window centered on the inner corner (1, 1)
        win = 3.5
        L = world_corridor_polygon(a).intersection(
            shp_box(1 - win, 1 - win, 1 + win, 1 + win))
        if L.geom_type == "Polygon":
            cx, cy = L.exterior.xy
            axA.fill(cx, cy, facecolor="#f1f1f1", edgecolor="#666666", linewidth=1.4)
        else:
            for g in L.geoms:
                cx, cy = g.exterior.xy
                axA.fill(cx, cy, facecolor="#f1f1f1", edgecolor="#666666", linewidth=1.4)
        axA.plot([1.0], [1.0], "o", color="#cc3300", markersize=4)
        axA.set_xlim(1 - win, 1 + win)
        axA.set_ylim(1 - win, 1 + win)
        _save(figA, f"aisle_{int(ad)}.png")

        # ----- sofa ------
        figS, axS = plt.subplots(figsize=(2.6, 2.6))
        axS.set_aspect("equal")
        axS.set_xticks([]); axS.set_yticks([])
        for s in axS.spines.values():
            s.set_visible(False)
        b = data[f"boundary_{i}"]
        if b.size:
            xs, ys = b[:, 0], b[:, 1]
            cx, cy = 0.5 * (xs.min() + xs.max()), 0.5 * (ys.min() + ys.max())
            sx, sy = xs - cx, ys - cy
            axS.add_patch(MplPolygon(np.column_stack([sx, sy]), closed=True,
                                     facecolor="#cfe2ff", edgecolor="#1f4fa8",
                                     linewidth=1.4))
        # Same scale across all sofas — use the largest's half-extent
        max_extent = max(xmax - xmin, ymax - ymin) / 2 + pad
        axS.set_xlim(-max_extent, max_extent)
        axS.set_ylim(-max_extent, max_extent)
        _save(figS, f"sofa_{int(ad)}.png")
        print(f"thumbs for alpha={ad}", flush=True)


def area_vs_angle(npz_path: str = "sweep_angles.npz") -> None:
    data = np.load(npz_path, allow_pickle=True)
    a = data["alphas_deg"].astype(float)
    A = data["areas"].astype(float)
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(a, A, "o-", color="#1f4fa8", linewidth=1.5, markersize=8,
            label="optimised sofa")
    # Reference curve: heuristic ~ C / (pi - alpha) scaled to match 90° value
    # (the agent's heuristic from the research note).
    aa = np.linspace(60, 175, 200)
    heuristic = 2.22 * (np.pi/2) / (np.pi - np.radians(aa))
    ax.plot(aa, heuristic, "--", color="#cc8800", linewidth=1.2,
            label=r"heuristic $\sim C / (\pi - \alpha)$")
    ax.axhline(np.pi / 2 + 2 / np.pi, color="#888888", linewidth=0.8,
               linestyle=":", label="Hammersley = π/2 + 2/π")
    ax.set_xlabel(r"bend angle $\alpha$ (degrees)")
    ax.set_ylabel(r"sofa area $S(\alpha)$")
    ax.set_title("Area of the optimised sofa vs. corridor bend angle")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_xlim(50, 175)
    ax.set_ylim(0, max(A) * 1.15)
    fig.tight_layout()
    _save(fig, "area_vs_angle.png")
    print("saved area_vs_angle.png")


def gerver_3d(npz_path: str = "sweep_angles.npz") -> None:
    """3D scene: extruded Gerver sofa sitting inside the L-tube.
    Square cross-section unit corridor (z in [0, 1]) bent 90 degrees.
    The extruded sofa is the trivial 3D lower bound, volume = 2.21953.
    """
    data = np.load(npz_path, allow_pickle=True)
    # index 2 is the 90-degree result
    b = data["boundary_2"]
    if b.size == 0:
        print("no 90 deg boundary; skipping 3D figure")
        return

    height = 1.0
    xs = b[:, 0].astype(float)
    ys = b[:, 1].astype(float)
    # Centre and reorient for a clean 3D look
    cx, cy = xs.mean(), ys.mean()
    xs -= cx; ys -= cy

    fig = plt.figure(figsize=(9.0, 5.5))
    ax = fig.add_subplot(111, projection="3d")

    # ---- L-tube (corridor outline) ----
    # Build a wireframe-ish L-tube along z, with floor at z = 0 and ceiling
    # at z = 1. Cross-section in (x, y) is the L-corridor; we draw two
    # outer L-rings at z = 0 and z = 1 plus a few verticals.
    arm_len = 4.5
    # We use the corridor in WORLD frame for the visualization (the
    # extruded sofa is at "midway" through the motion).
    L_xy = np.array([
        [-arm_len, 0], [1.0, 0],     # floor of arm 1 (in -x direction here)
        [1.0, -arm_len],             # outer wall going down
        [0.0, -arm_len], [0.0, 1.0], # ceiling of arm 2 (rotated)
        [-arm_len, 1.0],             # back to floor
        [-arm_len, 0],
    ], dtype=float)
    # NOTE: this draws a particular L-corridor; we will recentre below
    L_xy = L_xy - np.array([cx, cy])  # align with sofa centre
    bot_ring = np.column_stack([L_xy[:, 0], L_xy[:, 1], np.zeros(len(L_xy))])
    top_ring = np.column_stack([L_xy[:, 0], L_xy[:, 1], np.full(len(L_xy), height)])
    # Side walls of the L-tube (translucent)
    L_side_faces = []
    for i in range(len(L_xy) - 1):
        L_side_faces.append([bot_ring[i], bot_ring[i + 1],
                              top_ring[i + 1], top_ring[i]])
    ax.add_collection3d(Poly3DCollection(L_side_faces, alpha=0.18,
                                         facecolor="#cccccc", edgecolor="#666666",
                                         linewidth=0.6))
    # Add ground (floor) outline as a flat polygon for orientation
    ax.add_collection3d(Poly3DCollection(
        [list(bot_ring)], alpha=0.15,
        facecolor="#dddddd", edgecolor="#888888", linewidth=0.5))

    # ---- extruded Gerver sofa ----
    n = len(xs)
    bot = np.column_stack([xs, ys, np.zeros(n)])
    top = np.column_stack([xs, ys, np.full(n, height)])
    side_faces = []
    for i in range(n - 1):
        side_faces.append([bot[i], bot[i + 1], top[i + 1], top[i]])
    side_faces.append([bot[-1], bot[0], top[0], top[-1]])
    ax.add_collection3d(Poly3DCollection(side_faces, alpha=0.75,
                                         facecolor="#9ec8ff", edgecolor="#1f4fa8",
                                         linewidth=0.5))
    ax.add_collection3d(Poly3DCollection([list(bot)], alpha=0.75,
                                         facecolor="#7fb1ff", edgecolor="#1f4fa8",
                                         linewidth=0.5))
    ax.add_collection3d(Poly3DCollection([list(top)], alpha=0.75,
                                         facecolor="#bfd7ff", edgecolor="#1f4fa8",
                                         linewidth=0.5))

    ax.set_xlim(-3.5, 1.5); ax.set_ylim(-3.0, 2.0); ax.set_zlim(0, height + 0.4)
    ax.set_box_aspect((1.5, 1.5, 0.6))
    ax.view_init(elev=22, azim=-55)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title(
        f"Trivial 3D construction: extruded Gerver sofa inside the L-tube.\n"
        f"volume = 2D area × unit height = {float(data['areas'][2]):.4f} × 1 "
        f"= {float(data['areas'][2]):.4f}", fontsize=10)
    _save(fig, "gerver_3d.png")
    print("saved gerver_3d.png")


def tilted_3d_attempt(npz_path: str = "sweep_angles.npz") -> None:
    """Beyond the trivial extrusion: tilt the extruded Gerver sofa slightly
    about its long axis (the y-axis here). Demonstrates that tilt costs
    height clearance — a non-trivial 3D construction has to recover the
    lost volume by extending where the corridor still has room.
    """
    data = np.load(npz_path, allow_pickle=True)
    b = data["boundary_2"]
    if b.size == 0:
        return
    xs = b[:, 0].astype(float) - b[:, 0].mean()
    ys = b[:, 1].astype(float) - b[:, 1].mean()
    # tilt about y-axis: (x, z) rotate by tilt
    height = 1.0
    fig = plt.figure(figsize=(9.0, 4.0))
    for k, tilt_deg in enumerate([0, 10, 25]):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        t = math.radians(tilt_deg)
        n = len(xs)
        bot = np.array([[x, y, 0.0] for x, y in zip(xs, ys)])
        top = np.array([[x, y, height] for x, y in zip(xs, ys)])
        # Apply tilt around y-axis
        Ry = np.array([[math.cos(t), 0, math.sin(t)],
                       [0, 1, 0],
                       [-math.sin(t), 0, math.cos(t)]])
        bot = bot @ Ry.T
        top = top @ Ry.T
        # If any z > 1 (max corridor height), the sofa pokes through the
        # ceiling and would not fit
        clearance_top = top[:, 2].max()
        clearance_bot = bot[:, 2].min()
        valid = (clearance_top <= 1.0 + 1e-9) and (clearance_bot >= 0 - 1e-9)

        side_faces = []
        for i in range(n - 1):
            side_faces.append([bot[i], bot[i + 1], top[i + 1], top[i]])
        side_faces.append([bot[-1], bot[0], top[0], top[-1]])
        col = "#9ec8ff" if valid else "#ffb0b0"
        edgecol = "#1f4fa8" if valid else "#8b0000"
        ax.add_collection3d(Poly3DCollection(side_faces, alpha=0.7,
                                             facecolor=col, edgecolor=edgecol,
                                             linewidth=0.4))
        ax.add_collection3d(Poly3DCollection([list(bot)], alpha=0.7,
                                             facecolor=col, edgecolor=edgecol))
        ax.add_collection3d(Poly3DCollection([list(top)], alpha=0.7,
                                             facecolor=col, edgecolor=edgecol))
        ax.set_xlim(-2, 2); ax.set_ylim(-2.5, 2.5); ax.set_zlim(-0.5, 1.5)
        ax.set_box_aspect((1, 1.5, 0.6))
        ax.view_init(elev=18, azim=-50)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_title(f"tilt = {tilt_deg}°\n"
                     f"{'fits in the tube' if valid else 'pokes through ceiling/floor'}",
                     fontsize=9)
    fig.suptitle("Attempt: tilting the extruded Gerver sofa about its long axis",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "tilt_attempt.png")
    print("saved tilt_attempt.png")


def thumbs_for_table1(npz_path: str = "sweep_angles.npz") -> None:
    """Tiny shape thumbnails for Table 1."""
    # Unit square
    fig, ax = plt.subplots(figsize=(1.6, 1.6))
    ax.add_patch(MplPolygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True,
                            facecolor="#cfe2ff", edgecolor="#1f4fa8"))
    ax.set_aspect("equal")
    ax.set_xlim(-0.3, 1.3); ax.set_ylim(-0.3, 1.3)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    _save(fig, "thumb_unit_square.png")

    # Half disk
    fig, ax = plt.subplots(figsize=(1.6, 1.6))
    th = np.linspace(0, np.pi, 80)
    xs, ys = np.cos(th), np.sin(th)
    pts = np.column_stack([np.concatenate([[1], xs, [-1]]),
                           np.concatenate([[0], ys, [0]])])
    ax.add_patch(MplPolygon(pts, closed=True, facecolor="#cfe2ff", edgecolor="#1f4fa8"))
    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-0.3, 1.3)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    _save(fig, "thumb_half_disk.png")

    # Hammersley (analytic): rectangle + 2 quarter circles - semicircular notch
    r = 2.0 / np.pi
    fig, ax = plt.subplots(figsize=(2.6, 1.4))
    ax.set_aspect("equal")
    # outer boundary
    th_left = np.linspace(np.pi / 2, np.pi, 60)
    th_right = np.linspace(0, np.pi / 2, 60)
    outer = np.concatenate([
        np.column_stack([-r + np.cos(th_left), np.sin(th_left)]),
        np.column_stack([r + np.cos(th_right), np.sin(th_right)]),
    ])
    outer = np.vstack([
        [[r + 1, 0]], outer[::1],
        [[-r - 1, 0]], [[-r - 1, 0]], [[r + 1, 0]]
    ])
    # simpler: rectangle + arcs
    th = np.linspace(0, np.pi, 80)
    outer_top = []
    # Left quarter
    th_l = np.linspace(np.pi, np.pi / 2, 30)
    for t in th_l: outer_top.append([-r + np.cos(t), np.sin(t)])
    # Top straight (already covered by quarters meeting at y = 1)
    # Right quarter
    th_r = np.linspace(np.pi / 2, 0, 30)
    for t in th_r: outer_top.append([r + np.cos(t), np.sin(t)])
    bottom = [[r + 1, 0]]
    # notch (a semicircular hole, drawn as part of the boundary going inward)
    th_notch = np.linspace(0, np.pi, 60)
    notch = [[r * np.cos(t), r * np.sin(t)] for t in th_notch]
    boundary = outer_top + [bottom[0]] + notch[::-1] + [[-r - 1, 0]]
    boundary = np.array(boundary)
    ax.add_patch(MplPolygon(boundary, closed=True, facecolor="#cfe2ff", edgecolor="#1f4fa8"))
    ax.set_xlim(-2, 2); ax.set_ylim(-0.3, 1.3)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    _save(fig, "thumb_hammersley.png")

    # Gerver-class (our optimised) — copy the 90-degree boundary
    data = np.load(npz_path, allow_pickle=True)
    b = data["boundary_2"]
    if b.size:
        xs = b[:, 0] - b[:, 0].mean()
        ys = b[:, 1] - b[:, 1].mean()
        fig, ax = plt.subplots(figsize=(1.8, 2.6))
        ax.set_aspect("equal")
        ax.add_patch(MplPolygon(np.column_stack([xs, ys]), closed=True,
                                facecolor="#cfe2ff", edgecolor="#1f4fa8"))
        ax.set_xlim(xs.min() - 0.2, xs.max() + 0.2)
        ax.set_ylim(ys.min() - 0.2, ys.max() + 0.2)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        _save(fig, "thumb_gerver.png")
    print("table-1 thumbnails done")


if __name__ == "__main__":
    aisle_and_sofa_thumbs()
    area_vs_angle()
    gerver_3d()
    tilted_3d_attempt()
    thumbs_for_table1()
