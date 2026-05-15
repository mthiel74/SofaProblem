"""Dump everything the Wolfram figure-build script needs as plain CSV /
JSON.  The Wolfram script then has no dependency on Python beyond
having read these files once at the start of the post build."""

from __future__ import annotations

import json
import os
import numpy as np
from sofa_angle import (
    evaluate_trajectory_alpha,
    world_corridor_polygon,
    build_sofa_alpha,
)
from sofa_pipeline import (
    corridor_polygon as corridor_90,
    build_sofa as build_sofa_90,
    evaluate_trajectory as eval_traj_90,
)

OUT = "Wolfram/community/data"
os.makedirs(OUT, exist_ok=True)


def _save_pts(name, pts):
    np.savetxt(os.path.join(OUT, name), pts, delimiter=",", fmt="%.7f")


def export_2d_optimum():
    """The 90-degree optimised shape: sofa boundary, control points,
    and inner-corner trajectory densely sampled."""
    npz = np.load("sweep_angles.npz", allow_pickle=True)
    b = npz["boundary_2"].astype(float)
    ctrl = npz["controls_2"].astype(float)
    # densely-sampled corner trajectory
    eval_thetas = np.linspace(0.0, np.pi / 2.0, 240)
    corners = eval_traj_90(ctrl, eval_thetas)
    _save_pts("boundary_90.csv", b)
    _save_pts("controls_90.csv", ctrl)
    _save_pts("trajectory_90.csv", corners)
    np.savetxt(os.path.join(OUT, "thetas_90.csv"), eval_thetas,
               delimiter=",", fmt="%.7f")
    # Corridor polygon at the optimum at each theta, as a list of 6-vertex
    # polygons.  We'll stream them as a single file with "blank line"
    # separators.
    with open(os.path.join(OUT, "corridors_90.csv"), "w") as f:
        for theta, c in zip(eval_thetas, corners):
            poly = corridor_90(float(theta), float(c[0]), float(c[1]))
            for x, y in list(poly.exterior.coords)[:-1]:
                f.write(f"{x:.7f},{y:.7f}\n")
            f.write("\n")
    print("exported 2D optimum")


def export_angle_sweep():
    """Per-angle bend results: corridor outline, sofa boundary, trajectory."""
    data = np.load("sweep_angles.npz", allow_pickle=True)
    alphas = data["alphas_deg"].astype(int).tolist()
    areas = data["areas"].astype(float).tolist()
    summary = {"alphas_deg": alphas, "areas": areas}
    with open(os.path.join(OUT, "sweep_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    for i, ad in enumerate(alphas):
        # Corridor polygon (world frame) at this angle
        a = np.radians(float(ad))
        L = world_corridor_polygon(a)
        _save_pts(f"aisle_{int(ad)}.csv",
                   np.array(L.exterior.coords)[:-1])
        # Sofa boundary
        _save_pts(f"sofa_{int(ad)}.csv", data[f"boundary_{i}"].astype(float))
        # Control points
        _save_pts(f"controls_{int(ad)}.csv", data[f"controls_{i}"].astype(float))
        # Corner trajectory dense (1200 samples lives in eval_thetas / corners_dense)
        ev_th = data[f"eval_thetas_{i}"].astype(float)
        cd = data[f"corners_dense_{i}"].astype(float)
        np.savetxt(os.path.join(OUT, f"thetas_{int(ad)}.csv"), ev_th,
                   delimiter=",", fmt="%.7f")
        _save_pts(f"trajectory_{int(ad)}.csv", cd)
    print("exported angle sweep for", alphas)


if __name__ == "__main__":
    export_2d_optimum()
    export_angle_sweep()
    print(f"\nwrote {len(os.listdir(OUT))} files to {OUT}")
