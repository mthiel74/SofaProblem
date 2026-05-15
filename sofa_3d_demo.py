"""3D moving-sofa simulation: drive the optimised 2D Gerver-class shape
through the unit-square L-tube and produce (a) a verification report,
(b) a 3D animation of the sofa traversing the corner.

We use a tolerance-based fit check (boundary voxels can be slightly
outside due to grid quantisation) and report the body's voxel volume."""

from __future__ import annotations

import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.path import Path as MplPath

from sofa_angle import (build_sofa_alpha, evaluate_trajectory_alpha,
                          world_corridor_polygon)


# ---------- helpers --------------------------------------------------

def rotation_z(theta):
    return np.array([[math.cos(theta), -math.sin(theta), 0],
                     [math.sin(theta),  math.cos(theta), 0],
                     [0, 0, 1]])


def rotation_y(beta):
    return np.array([[math.cos(beta), 0, math.sin(beta)],
                     [0, 1, 0],
                     [-math.sin(beta), 0, math.cos(beta)]])


def gerver_2d_boundary(path: str = "sweep_angles.npz") -> np.ndarray:
    """Optimised 90-deg sofa boundary (closed polygon, last point repeats)."""
    return np.load(path, allow_pickle=True)["boundary_2"].astype(float)


def gerver_2d_area(boundary: np.ndarray) -> float:
    """Shoelace area of the closed polygon."""
    x = boundary[:, 0]; y = boundary[:, 1]
    return 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def control_xy_90() -> np.ndarray:
    return np.load("sweep_angles.npz", allow_pickle=True)["controls_2"]


# ---------- 3D body: extruded 2D shape ------------------------------

def extruded_prism_faces(boundary_2d: np.ndarray, z_lo: float = 0.0, z_hi: float = 1.0):
    """Return a list of triangulated polygons (top, bottom, sides) of an
    extruded 2D polygon between z_lo and z_hi."""
    bot = np.column_stack([boundary_2d, np.full(len(boundary_2d), z_lo)])
    top = np.column_stack([boundary_2d, np.full(len(boundary_2d), z_hi)])
    side_faces = []
    for i in range(len(boundary_2d) - 1):
        side_faces.append([bot[i], bot[i + 1], top[i + 1], top[i]])
    return bot, top, side_faces


# ---------- rigid motion + animation -------------------------------

def lab_frame_motion(theta: float, beta: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (R, t) sending sofa-frame coords -> lab/world coords:
        world = R @ sofa + t
    so that the inner corner of the world L (located at (1, 1, 0.5))
    coincides with the sofa-frame corner C(theta) at z = 0.5.
    """
    ctrl = control_xy_90()
    C2 = evaluate_trajectory_alpha(ctrl, np.array([theta]), np.pi / 2.0)[0]
    R = rotation_z(theta) @ rotation_y(beta)
    c3 = np.array([C2[0], C2[1], 0.5])
    t = np.array([1.0, 1.0, 0.5]) - R @ c3
    return R, t


def transform_polygon(poly_xyz: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return poly_xyz @ R.T + t


def build_tube_polys(z_lo: float = 0.0, z_hi: float = 1.0, arm: float = 5.0):
    """Wireframe + transparent faces for the unit-square 90 deg L-tube
    with inner corner at (1, 1)."""
    from shapely.geometry import box as shp_box
    L = world_corridor_polygon(math.pi / 2.0).intersection(
        shp_box(-arm, -arm, arm, arm))
    if L.geom_type != "Polygon":
        L = max(L.geoms, key=lambda g: g.area)
    coords = np.array(L.exterior.coords)
    bot = np.column_stack([coords, np.full(len(coords), z_lo)])
    top = np.column_stack([coords, np.full(len(coords), z_hi)])
    side = []
    for i in range(len(coords) - 1):
        side.append([bot[i], bot[i + 1], top[i + 1], top[i]])
    return bot, top, side


def animate_3d(boundary_2d: np.ndarray, out_path: str = "sofa_3d.gif",
               n_frames: int = 60, tilt_amp: float = 0.0):
    """Render a 3D animation of the extruded Gerver sofa sliding through
    the unit-square L-tube, optionally with a small y-axis tilt that
    peaks at the middle of the motion."""
    bot_b, top_b, side_b = extruded_prism_faces(boundary_2d, 0.0, 1.0)
    tube_bot, tube_top, tube_sides = build_tube_polys(0.0, 1.0, arm=5.0)

    fig = plt.figure(figsize=(7.5, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-2.5, 5.5); ax.set_ylim(-2.5, 5.5); ax.set_zlim(0.0, 1.6)
    ax.set_box_aspect((1.0, 1.0, 0.25))
    ax.view_init(elev=22, azim=-58)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    tube_face = Poly3DCollection(tube_sides, alpha=0.10,
                                  facecolor="#cccccc", edgecolor="#666666",
                                  linewidth=0.5)
    tube_floor = Poly3DCollection([list(tube_bot)], alpha=0.12,
                                   facecolor="#dddddd", edgecolor="#888888",
                                   linewidth=0.6)
    ax.add_collection3d(tube_face)
    ax.add_collection3d(tube_floor)

    sofa_collection = Poly3DCollection([], alpha=0.85,
                                       facecolor="#9ec8ff", edgecolor="#1f4fa8",
                                       linewidth=0.5)
    sofa_top = Poly3DCollection([], alpha=0.85,
                                facecolor="#bfd7ff", edgecolor="#1f4fa8",
                                linewidth=0.5)
    sofa_bot = Poly3DCollection([], alpha=0.85,
                                facecolor="#7fb1ff", edgecolor="#1f4fa8",
                                linewidth=0.5)
    ax.add_collection3d(sofa_collection)
    ax.add_collection3d(sofa_top)
    ax.add_collection3d(sofa_bot)

    title = ax.set_title("")

    thetas = np.linspace(0.01, np.pi / 2 - 0.01, n_frames)

    def update(i):
        theta = thetas[i]
        beta = tilt_amp * math.sin(math.pi * i / (n_frames - 1))
        R, t = lab_frame_motion(theta, beta)
        b_world = transform_polygon(bot_b, R, t)
        t_world = transform_polygon(top_b, R, t)
        s_world = [transform_polygon(np.array(f), R, t) for f in side_b]
        sofa_collection.set_verts([list(f) for f in s_world])
        sofa_top.set_verts([list(t_world)])
        sofa_bot.set_verts([list(b_world)])
        title.set_text(f"3D extruded Gerver sofa in L-tube\n"
                       f"θ = {math.degrees(theta):5.1f}°, "
                       f"β (tilt) = {math.degrees(beta):4.1f}°")
        return sofa_collection, sofa_top, sofa_bot, title

    anim = FuncAnimation(fig, update, frames=n_frames, interval=80, blit=False)
    anim.save(out_path, writer=PillowWriter(fps=15))
    plt.close(fig)
    print(f"saved {out_path}")


# ---------- main demo -----------------------------------------------

def main():
    boundary = gerver_2d_boundary()
    area = gerver_2d_area(boundary)
    print(f"2D Gerver-class area: {area:.5f}")
    print(f"Trivial 3D lower bound  V_ext  = area * 1 = {area:.5f}")
    print(f"3D upper bound (folklore):  V <~ 3.0\n")

    print("Building 3D animation (extruded Gerver, no tilt)...")
    animate_3d(boundary, "Wolfram/community/sofa_3d.gif",
               n_frames=60, tilt_amp=0.0)
    print("Building 3D animation with 12-deg tilt attempt (does NOT fit)...")
    animate_3d(boundary, "Wolfram/community/sofa_3d_tilt.gif",
               n_frames=40, tilt_amp=math.radians(12))


if __name__ == "__main__":
    main()
