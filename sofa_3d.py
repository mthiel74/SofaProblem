"""A first numerical attack on the 3D moving-sofa problem.

The 3D L-tube is the Cartesian product (2D L-corridor) x [0, 1] in z.
A candidate 3D sofa is a rigid body S subset of R^3; admissibility
asks for a continuous SE(3) motion g(t) keeping g(t) S contained in
the tube. We check candidates numerically by:

  1. voxelising the tube and the candidate body,
  2. discretising the motion as a sequence of (R, t),
  3. for each (R, t), verifying that the rotated/translated voxels
     of the body are a subset of the tube voxels,
  4. reporting volume = (#body voxels) * (voxel size)^3.

The trivial extrusion S_2D x [0,1] should give V = area(S_2D) = 2.21953.
We then try a few non-trivial constructions (rigid tilt; modified
prism with a "fin" added in a 3D-accessible region) and report.

This is intentionally a coarse numerical probe of an open problem;
the framework is what matters more than the precise voxel resolution.
"""

from __future__ import annotations

import math
import numpy as np
from sofa_angle import build_sofa_alpha, evaluate_trajectory_alpha


# --------------------- 2D Gerver-class sofa (cached) ---------------------

def load_2d_gerver_boundary(path: str = "sweep_angles.npz") -> np.ndarray:
    """Return the boundary polygon of our optimised 90-degree sofa."""
    data = np.load(path, allow_pickle=True)
    return data["boundary_2"].astype(float)


def gerver_2d_membership(xy_grid: np.ndarray, boundary: np.ndarray) -> np.ndarray:
    """Boolean mask, True where (x, y) is inside the 2D Gerver polygon."""
    from matplotlib.path import Path
    p = Path(boundary)
    return p.contains_points(xy_grid).astype(bool)


# --------------------- Voxel grid ---------------------------------------

class VoxelGrid:
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, res):
        self.res = float(res)
        self.xmin = xmin; self.xmax = xmax
        self.ymin = ymin; self.ymax = ymax
        self.zmin = zmin; self.zmax = zmax
        self.nx = int(math.ceil((xmax - xmin) / res))
        self.ny = int(math.ceil((ymax - ymin) / res))
        self.nz = int(math.ceil((zmax - zmin) / res))

    def cell_centers(self):
        xs = self.xmin + (np.arange(self.nx) + 0.5) * self.res
        ys = self.ymin + (np.arange(self.ny) + 0.5) * self.res
        zs = self.zmin + (np.arange(self.nz) + 0.5) * self.res
        return xs, ys, zs

    def info(self):
        return f"VoxelGrid {self.nx}x{self.ny}x{self.nz}, res={self.res}"


# --------------------- 3D corridor as a voxel mask ------------------------

def build_tube_mask(grid: VoxelGrid, alpha_deg: float = 90.0) -> np.ndarray:
    """L-tube (2D L-corridor) x [0, 1]. Inner corner at (1, 1) in (x, y);
    z floor at 0 and ceiling at 1.  Returns a (nx, ny, nz) bool array."""
    from sofa_angle import world_corridor_polygon
    from matplotlib.path import Path
    a = math.radians(alpha_deg)
    poly = world_corridor_polygon(a)
    xs, ys, zs = grid.cell_centers()
    # 2D mask first
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([XX.ravel(), YY.ravel()])
    # Clip the corridor to a tractable region — use the bounding box of
    # the voxel grid.
    from shapely.geometry import box as shp_box
    clip = poly.intersection(shp_box(grid.xmin, grid.ymin, grid.xmax, grid.ymax))
    if clip.is_empty:
        return np.zeros((grid.nx, grid.ny, grid.nz), dtype=bool)
    if clip.geom_type == "Polygon":
        m = Path(np.array(clip.exterior.coords)).contains_points(pts)
    else:
        # union of polygons -> OR of memberships
        m = np.zeros(pts.shape[0], dtype=bool)
        for g in clip.geoms:
            m = m | Path(np.array(g.exterior.coords)).contains_points(pts)
    mask2d = m.reshape(grid.nx, grid.ny)
    z_in = (zs >= 0.0) & (zs <= 1.0)
    return mask2d[:, :, None] & z_in[None, None, :]


# --------------------- Candidate body voxels ----------------------------

def trivial_extrusion_body(grid: VoxelGrid, boundary_2d: np.ndarray,
                            z_max: float = 1.0) -> np.ndarray:
    """The 2D Gerver shape extruded along z from 0 to z_max."""
    from matplotlib.path import Path
    xs, ys, zs = grid.cell_centers()
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([XX.ravel(), YY.ravel()])
    in_2d = Path(boundary_2d).contains_points(pts).reshape(grid.nx, grid.ny)
    z_in = (zs >= 0.0) & (zs <= z_max)
    return in_2d[:, :, None] & z_in[None, None, :]


# --------------------- Rigid motion ------------------------------------

def rotation_3d(theta: float, beta: float = 0.0) -> np.ndarray:
    """Rotation around the (vertical) z-axis by theta, then around the
    long axis (here taken as the y-axis) by beta. The y-axis tilt is
    the extra 3D degree of freedom we're exercising."""
    Rz = np.array([[math.cos(theta), -math.sin(theta), 0],
                   [math.sin(theta),  math.cos(theta), 0],
                   [0, 0, 1]])
    Ry = np.array([[math.cos(beta), 0, math.sin(beta)],
                   [0, 1, 0],
                   [-math.sin(beta), 0, math.cos(beta)]])
    return Rz @ Ry


def check_body_in_tube(body_voxel_coords: np.ndarray,
                        R: np.ndarray, t: np.ndarray,
                        tube_mask: np.ndarray,
                        grid: VoxelGrid) -> tuple[bool, int]:
    """Apply rigid motion (R, t) to body voxel CENTRE coordinates; check
    that every transformed centre falls in a tube voxel.

    Returns (all_inside, n_outside).
    """
    moved = body_voxel_coords @ R.T + t
    ix = np.floor((moved[:, 0] - grid.xmin) / grid.res).astype(int)
    iy = np.floor((moved[:, 1] - grid.ymin) / grid.res).astype(int)
    iz = np.floor((moved[:, 2] - grid.zmin) / grid.res).astype(int)
    in_bounds = (ix >= 0) & (ix < grid.nx) & \
                (iy >= 0) & (iy < grid.ny) & \
                (iz >= 0) & (iz < grid.nz)
    inside_tube = np.zeros(len(moved), dtype=bool)
    inside_tube[in_bounds] = tube_mask[ix[in_bounds], iy[in_bounds], iz[in_bounds]]
    n_out = int((~inside_tube).sum())
    return n_out == 0, n_out


# --------------------- Pipeline ----------------------------------------

def sweep_body_through_corner(body_voxel_coords: np.ndarray,
                              motion: list[tuple[np.ndarray, np.ndarray]],
                              grid: VoxelGrid, tube_mask: np.ndarray,
                              verbose: bool = False) -> dict:
    """Apply each motion step to the body voxels; report worst-case
    out-of-tube count."""
    worst = 0
    fits_every_step = True
    for i, (R, t) in enumerate(motion):
        ok, n_out = check_body_in_tube(body_voxel_coords, R, t,
                                       tube_mask, grid)
        if n_out > worst:
            worst = n_out
        if not ok:
            fits_every_step = False
            if verbose:
                print(f"  step {i}: {n_out} voxels outside tube")
                break
    return {"fits": fits_every_step, "worst_n_outside": worst}


def make_corner_motion(n_steps: int = 40, tilt_amp: float = 0.0) -> list:
    """Generate a 2D corner motion derived from the 90-degree optimised
    trajectory, with an optional small y-axis tilt that oscillates
    through the corner (0 at the ends, ±tilt_amp in the middle)."""
    data = np.load("sweep_angles.npz", allow_pickle=True)
    ctrl = data["controls_2"]  # 90 deg controls
    theta_max = np.pi / 2
    thetas = np.linspace(0.0, theta_max, n_steps)
    corners = evaluate_trajectory_alpha(ctrl, thetas, theta_max)
    motions = []
    for i, (theta, c) in enumerate(zip(thetas, corners)):
        # The corner moves over time; in sofa-fixed frame the corridor
        # rotates by -theta around the world inner corner (1, 1). Equivalent
        # rigid motion of the sofa in the world frame:
        # R(theta) about world's z-axis, translate so that the sofa's
        # frame origin (which sits at C(theta) in sofa coords) lands on the
        # world inner corner (1, 1).
        # In the sofa-fixed frame: g_sofa^(-1) takes world to sofa.
        # In the lab frame this is the rotation that places the sofa.
        beta = tilt_amp * math.sin(math.pi * i / (n_steps - 1))
        R = rotation_3d(theta, beta)
        # Translate so that the sofa-frame corner c lands at world (1, 1, 0.5)
        # (the corridor centre). Apply R^T to map sofa-frame to world.
        c3 = np.array([c[0], c[1], 0.5])
        t = np.array([1.0, 1.0, 0.5]) - R @ c3
        motions.append((R, t))
    return motions


# --------------------- Volume scoring -----------------------------------

def body_volume(body_mask: np.ndarray, grid: VoxelGrid) -> float:
    return float(body_mask.sum()) * grid.res ** 3


def body_voxel_centers(body_mask: np.ndarray, grid: VoxelGrid) -> np.ndarray:
    xs, ys, zs = grid.cell_centers()
    ix, iy, iz = np.where(body_mask)
    return np.column_stack([xs[ix], ys[iy], zs[iz]])


# --------------------- Demonstrations -----------------------------------

def main():
    boundary = load_2d_gerver_boundary()
    print(f"2D Gerver boundary: {boundary.shape[0]} vertices")

    # Build a voxel grid spanning the corridor neighbourhood
    res = 0.06
    grid = VoxelGrid(xmin=-4.0, xmax=4.0, ymin=-4.0, ymax=4.0,
                     zmin=-0.2, zmax=1.2, res=res)
    print(grid.info())

    tube = build_tube_mask(grid)
    print(f"tube voxels: {int(tube.sum())} of {tube.size}")

    # Candidate 1: trivial extrusion
    body = trivial_extrusion_body(grid, boundary, z_max=1.0)
    centers = body_voxel_centers(body, grid)
    V = body_volume(body, grid)
    print(f"\nCandidate 1: trivial Gerver extrusion (height 1)")
    print(f"  voxel volume = {V:.4f}  (target: 2D area = 2.21968)")

    # Sweep through the corner with 2D motion only
    motion = make_corner_motion(n_steps=30, tilt_amp=0.0)
    res_pure = sweep_body_through_corner(centers, motion, grid, tube)
    print(f"  pure 2D motion:  fits={res_pure['fits']}, "
          f"worst_n_outside={res_pure['worst_n_outside']}")

    # Same body with a small y-axis tilt over the motion
    for tilt_deg in (5, 10, 20):
        m = make_corner_motion(n_steps=30,
                                tilt_amp=math.radians(tilt_deg))
        r = sweep_body_through_corner(centers, m, grid, tube)
        print(f"  with {tilt_deg} deg tilt: fits={r['fits']}, "
              f"worst_n_outside={r['worst_n_outside']}")

    # Candidate 2: extruded Gerver with a small wedge added on top
    # (a "horn" extending above z = 1 — should violate the ceiling
    # under pure 2D motion, but might survive under a tilt).
    horn = body.copy()
    xs, ys, zs = grid.cell_centers()
    # Add a small slab at z in [1, 1.1] for the inner portion of the Gerver
    # base.  Use the same 2D mask but trimmed by x.
    from matplotlib.path import Path
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([XX.ravel(), YY.ravel()])
    in2d = Path(boundary).contains_points(pts).reshape(grid.nx, grid.ny)
    xs_2d = XX[:, 0]
    xmask = xs_2d >= 0.5  # only points with x >= 0.5 in body coords
    extra = in2d & xmask[:, None]
    z_horn = (zs >= 1.0) & (zs <= 1.1)
    horn = body | (extra[:, :, None] & z_horn[None, None, :])
    centers_h = body_voxel_centers(horn, grid)
    Vh = body_volume(horn, grid)
    print(f"\nCandidate 2: Gerver prism + thin 'horn' on top (x>0.5, z in [1, 1.1])")
    print(f"  voxel volume = {Vh:.4f}")
    for tilt_deg in (0, 5, 10):
        m = make_corner_motion(n_steps=30,
                                tilt_amp=math.radians(tilt_deg))
        r = sweep_body_through_corner(centers_h, m, grid, tube)
        print(f"  tilt={tilt_deg} deg: fits={r['fits']}, "
              f"worst_n_outside={r['worst_n_outside']}")


if __name__ == "__main__":
    main()
