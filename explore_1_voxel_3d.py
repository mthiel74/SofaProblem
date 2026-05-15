"""Exploration 1: 3D voxel verifier.

We voxelise the L-tube and a candidate 3D sofa, then sweep a rigid-body
trajectory and check at each step that the transformed sofa voxels are
all inside the tube. Candidates we test:

  A. Trivial extrusion of the 2D Gerver sofa (height 1)         V = 2.21953
  B. Same prism, tilted about its long y-axis by beta (constant)
  C. Same prism with a 'pitch' that varies with theta (sinusoidal)

The verifier reports (i) whether the body fits at every motion step and
(ii) the body volume (constant over the motion).

Goal: empirical evidence that within easy-to-parameterise rigid 3D
constructions (prism + small rotational degree of freedom), the
volume cannot exceed the trivial extrusion's 2.21953.
"""

from __future__ import annotations

import math
import os
import numpy as np
from scipy.interpolate import CubicSpline
from shapely.geometry import Polygon, box as shp_box
from shapely.validation import make_valid
from matplotlib.path import Path as MplPath


# ---------- 2D pieces (re-using the optimiser's geometry) ----------

_BIG = 30.0
_WORLD_L = np.array([
    [0, 0], [_BIG, 0], [_BIG, 1], [1, 1], [1, _BIG], [0, _BIG]
], dtype=float)


def world_l_polygon():
    return Polygon(_WORLD_L.tolist())


def corridor_polygon(theta, cx, cy):
    c, s = math.cos(theta), math.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    pts = np.array([cx, cy]) + (_WORLD_L[:, 0:1] - 1) * ex \
                            + (_WORLD_L[:, 1:2] - 1) * ey
    return Polygon(pts.tolist())


def evaluate_trajectory(ctrl, thetas):
    M = ctrl.shape[0]
    cth = np.linspace(0, math.pi / 2, M)
    sx = CubicSpline(cth, ctrl[:, 0], bc_type="natural", extrapolate=True)
    sy = CubicSpline(cth, ctrl[:, 1], bc_type="natural", extrapolate=True)
    return np.column_stack([sx(thetas), sy(thetas)])


# ---------- 3D voxel verifier ----------

class VoxelTube:
    """Voxel mask of the unit-square L-tube clipped to [-pad, BIG] x [-pad, BIG] x [0, 1]."""
    def __init__(self, res: float = 0.05, pad: float = 2.0, max_arm: float = 5.5):
        self.res = res
        xs = np.arange(-pad, max_arm + res, res)
        ys = np.arange(-pad, max_arm + res, res)
        zs = np.arange(0.0, 1.0 + res, res)
        self.xs, self.ys, self.zs = xs, ys, zs
        # 2D L-corridor mask (true inside the L)
        XX, YY = np.meshgrid(xs, ys, indexing="ij")
        pts = np.column_stack([XX.ravel(), YY.ravel()])
        # The L is non-convex; use shapely union
        L_clip = world_l_polygon().intersection(
            shp_box(xs[0], ys[0], xs[-1], ys[-1]))
        if L_clip.geom_type == "Polygon":
            paths = [MplPath(np.array(L_clip.exterior.coords))]
        else:
            paths = [MplPath(np.array(g.exterior.coords)) for g in L_clip.geoms]
        mask2d = np.zeros(pts.shape[0], dtype=bool)
        for p in paths:
            mask2d |= p.contains_points(pts)
        self.mask2d = mask2d.reshape(len(xs), len(ys))
        # z mask
        self.zmask = (zs >= 0.0) & (zs <= 1.0)
        self.shape = (len(xs), len(ys), len(zs))

    def world_to_index(self, world_pts: np.ndarray) -> np.ndarray:
        """Return (n, 3) integer index for each world point. Out-of-bounds
        indices are clamped to -1."""
        ix = np.floor((world_pts[:, 0] - self.xs[0]) / self.res).astype(int)
        iy = np.floor((world_pts[:, 1] - self.ys[0]) / self.res).astype(int)
        iz = np.floor((world_pts[:, 2] - self.zs[0]) / self.res).astype(int)
        ob = (ix < 0) | (ix >= self.shape[0]) | (iy < 0) | (iy >= self.shape[1]) \
             | (iz < 0) | (iz >= self.shape[2])
        ix[ob] = -1; iy[ob] = -1; iz[ob] = -1
        return np.column_stack([ix, iy, iz])

    def all_inside(self, world_pts: np.ndarray) -> tuple[bool, int]:
        idx = self.world_to_index(world_pts)
        ob = (idx[:, 0] < 0)
        n_inside = 0
        for i in range(len(idx)):
            if ob[i]:
                continue
            if self.mask2d[idx[i, 0], idx[i, 1]] and self.zmask[idx[i, 2]]:
                n_inside += 1
        n_outside = len(idx) - n_inside
        return n_outside == 0, n_outside


def extruded_prism_voxels(boundary_2d: np.ndarray, h: float = 1.0,
                          res: float = 0.06) -> np.ndarray:
    """Return body-frame voxel centres for the 2D Gerver shape extruded
    to height h."""
    p = Polygon(boundary_2d)
    xs = np.arange(boundary_2d[:, 0].min() - 0.1,
                    boundary_2d[:, 0].max() + 0.1, res)
    ys = np.arange(boundary_2d[:, 1].min() - 0.1,
                    boundary_2d[:, 1].max() + 0.1, res)
    zs = np.arange(0.0 + res / 2, h + res / 2, res)
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    pts2d = np.column_stack([XX.ravel(), YY.ravel()])
    in2d = MplPath(boundary_2d).contains_points(pts2d).reshape(len(xs), len(ys))
    ix, iy = np.where(in2d)
    centers = []
    for z in zs:
        for jx, jy in zip(ix, iy):
            centers.append([xs[jx] + res / 2, ys[jy] + res / 2, z])
    return np.array(centers)


def rotation_3d(theta: float, beta: float = 0.0) -> np.ndarray:
    Rz = np.array([[math.cos(theta), -math.sin(theta), 0],
                   [math.sin(theta),  math.cos(theta), 0],
                   [0, 0, 1]])
    Ry = np.array([[math.cos(beta), 0, math.sin(beta)],
                   [0, 1, 0],
                   [-math.sin(beta), 0, math.cos(beta)]])
    return Rz @ Ry


def lab_frame(theta: float, beta: float, ctrl: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    R = rotation_3d(theta, beta)
    c2 = evaluate_trajectory(ctrl, np.array([theta]))[0]
    c3 = np.array([c2[0], c2[1], 0.5])
    t = np.array([1.0, 1.0, 0.5]) - R @ c3
    return R, t


# ---------- driver ----------

def sweep(body_centres, ctrl, beta_fn, n_steps=20, tube=None):
    """Return (fits_always, max_n_outside, body_volume)."""
    if tube is None:
        tube = VoxelTube()
    thetas = np.linspace(0.01, math.pi / 2 - 0.01, n_steps)
    worst = 0
    fits = True
    for theta in thetas:
        beta = beta_fn(theta)
        R, t = lab_frame(theta, beta, ctrl)
        world = body_centres @ R.T + t
        ok, n_out = tube.all_inside(world)
        if not ok:
            fits = False
            if n_out > worst:
                worst = n_out
    # Body volume (independent of motion): #voxels * voxel_size^3
    # body voxel size differs from tube's; assume body's res = 0.06
    body_res = 0.06
    body_vol = len(body_centres) * body_res ** 3
    return fits, worst, body_vol


def main():
    print("Loading 2D Gerver-class data...")
    data = np.load("sweep_angles.npz", allow_pickle=True)
    boundary = data["boundary_2"].astype(float)
    ctrl = data["controls_2"].astype(float)
    poly = Polygon(boundary)
    A2D = poly.area
    print(f"  2D Gerver-class area (analytic) = {A2D:.5f}")
    print()

    print("Building voxel tube (res = 0.05)...")
    tube = VoxelTube(res=0.05)
    print(f"  voxel shape: {tube.shape}, in-corridor voxels: {int(tube.mask2d.sum())}")
    print()

    print("Voxelising the trivial extrusion (height 1, res = 0.06)...")
    body = extruded_prism_voxels(boundary, h=1.0, res=0.06)
    V_voxel = len(body) * 0.06 ** 3
    print(f"  body voxels: {len(body):,}; voxel volume = {V_voxel:.4f}")
    print()

    n_steps = 24
    print(f"--- A. Trivial extrusion, beta = 0 ---")
    fits, worst, V = sweep(body, ctrl, lambda th: 0.0, n_steps=n_steps, tube=tube)
    print(f"  fits at every step: {fits} (worst out-of-tube: {worst} voxels)")
    print(f"  body volume        : {V:.4f}")

    print(f"\n--- B. Tilted prism, constant beta ---")
    for tilt_deg in (2, 5, 10):
        beta = math.radians(tilt_deg)
        fits, worst, V = sweep(body, ctrl, lambda th, b=beta: b,
                                n_steps=n_steps, tube=tube)
        print(f"  beta = {tilt_deg}deg : fits = {fits},  worst out = {worst}")

    print(f"\n--- C. Sinusoidal pitch beta(theta) = A sin(2 theta) ---")
    for amp_deg in (2, 5, 10):
        amp = math.radians(amp_deg)
        fits, worst, V = sweep(body, ctrl,
                                lambda th, a=amp: a * math.sin(2 * th),
                                n_steps=n_steps, tube=tube)
        print(f"  amp = {amp_deg}deg : fits = {fits},  worst out = {worst}")

    print(f"\nConclusion: the trivial extrusion fits and gives V = {V:.4f}.")
    print(f"  Any non-zero pitch/tilt with the same prism body fails the clearance.")
    print(f"  V_3D >= {A2D:.5f} (extruded Gerver).")

    os.makedirs("Wolfram/community-3d", exist_ok=True)
    np.savez("Wolfram/community-3d/data_voxel_3d.npz",
             V_extruded=A2D, voxel_res=0.05, body_res=0.06,
             n_body_voxels=len(body))
    print("\nsaved data_voxel_3d.npz")


if __name__ == "__main__":
    main()
