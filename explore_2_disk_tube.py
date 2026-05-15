"""Exploration 2: 3D L-tube with a DISK cross-section.

Unit-diameter circular cross-section instead of unit square. At height
z in [-1/2, 1/2] above the centre, the available width perpendicular to
the motion direction is
    w(z) = 2 sqrt(1/4 - z^2).

For each z-slice of a candidate 3D sofa, that slice has to fit in an
L-corridor of width w(z). The largest 2D sofa for L-width w is A(w),
which we tabulated in explore_4_variable_width.

Two bounds:

  V_prism :   prismatic sofa (constant cross-section over z in [-h/2, h/2]).
              The body width is the narrowest available, so
                  V_prism = max_h A(sqrt(1 - h^2)) * h.

  V_slices :  slice-by-slice upper bound (each slice independently optimised,
              ignoring rigid-body constraints):
                  V_slices = integral over z in [-1/2, 1/2] of A(w(z)) dz
                           = 2 * integral over z in [0, 1/2] of A(2 sqrt(1/4 - z^2)) dz.
              This is NOT achievable by a rigid body in general (slices
              have different optimal motions), but bounds the disk-tube
              sofa volume from above.

We report both numbers and the v_max - v_min gap, which quantifies how
much the rigid-body constraint costs in the disk tube.
"""

from __future__ import annotations

import math
import os
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator


def main():
    data = np.load("Wolfram/community-3d/data_variable_width.npz")
    widths = data["widths"].astype(float)
    areas  = data["areas"].astype(float)
    # Cubic monotone interpolation; A(0) extrapolated to 0.
    widths_ext = np.concatenate([[0.0], widths])
    areas_ext  = np.concatenate([[0.0], areas])
    order = np.argsort(widths_ext)
    A = PchipInterpolator(widths_ext[order], areas_ext[order],
                          extrapolate=False)

    print("variable-width data (from exploration 4):")
    for w, a in sorted(zip(widths, areas)):
        print(f"  A({w:.2f}) = {a:.5f}")
    print()

    # --- V_prism: prismatic sofa volume ---
    def prism_vol(h):
        if h <= 0 or h >= 1: return 0.0
        w = math.sqrt(1 - h * h)
        if w > widths.max(): w = widths.max()
        return float(A(w)) * h

    hs = np.linspace(0.0, 1.0, 1001)
    vs = [prism_vol(h) for h in hs]
    h_star = hs[int(np.argmax(vs))]
    V_prism = max(vs)
    w_star = math.sqrt(1 - h_star * h_star)
    print(f"V_prism (best constant cross-section):")
    print(f"  h*    = {h_star:.4f}")
    print(f"  w*    = {w_star:.4f}")
    print(f"  A(w*) = {float(A(min(w_star, widths.max()))):.4f}")
    print(f"  V_prism = h * A(w) = {V_prism:.4f}")
    print()

    # --- V_slices: slice-by-slice upper bound ---
    def w_of_z(z):
        if z**2 >= 0.25:  return 0.0
        return 2.0 * math.sqrt(0.25 - z * z)
    def slice_area(z):
        w = w_of_z(z)
        if w <= 0: return 0.0
        if w > widths.max(): w = widths.max()
        return float(A(w))
    V_slices, err = quad(slice_area, -0.5, 0.5, limit=200)
    print(f"V_slices (slice-by-slice upper bound, no rigid-body constraint):")
    print(f"  V_slices = integral of A(w(z)) dz = {V_slices:.4f}  (err {err:.2e})")
    print()

    # --- comparison ---
    print("disk-tube 3D bounds:")
    print(f"  V_prism   = {V_prism:.4f}  (achievable rigid-body lower bound)")
    print(f"  V_slices  = {V_slices:.4f}  (slice-by-slice upper bound)")
    print(f"  ratio     = {V_prism / V_slices:.4f}")
    print(f"  for comparison: square-tube trivial extrusion = 2.21953")

    np.savez("Wolfram/community-3d/data_disk_tube.npz",
             hs=hs, prism_vols=np.array(vs),
             h_star=h_star, V_prism=V_prism, w_star=w_star,
             V_slices=V_slices,
             widths_data=widths, areas_data=areas)
    print("\nsaved data_disk_tube.npz")


if __name__ == "__main__":
    main()
