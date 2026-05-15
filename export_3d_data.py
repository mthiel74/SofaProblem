"""Export numerical results from the 5 explorations as CSV for the
Wolfram figure builder."""

import os
import numpy as np

OUT = "Wolfram/community-3d"
os.makedirs(OUT, exist_ok=True)


def save(name, arr):
    np.savetxt(os.path.join(OUT, name), arr, delimiter=",", fmt="%.6f")


def main():
    # ---- #4 variable width ----
    d = np.load(f"{OUT}/data_variable_width.npz")
    widths = d["widths"]; areas = d["areas"]
    save("data_widths.csv", widths[:, None])
    save("data_areas.csv", areas[:, None])
    # individual width boundaries for the gallery
    for i, w in enumerate(widths):
        save(f"data_width_{str(float(w)).replace('.', 'p')}.csv",
              d[f"boundary_{i}"].astype(float))

    # ---- #3 S-corridor ----
    d = np.load(f"{OUT}/data_s_corridor.npz")
    save("data_s_boundary.csv", d["boundary"].astype(float))

    # ---- #5 convex ----
    d = np.load(f"{OUT}/data_convex.npz")
    save("data_convex_boundary.csv", d["boundary"].astype(float))
    # Also save Gerver-class boundary for the side-by-side comparison
    g = np.load("sweep_angles.npz", allow_pickle=True)
    save("data_gerver_boundary.csv", g["boundary_2"].astype(float))

    # ---- #2 disk tube ----
    d = np.load(f"{OUT}/data_disk_tube.npz")
    save("data_disk_hs.csv", d["hs"][:, None])
    save("data_disk_vs.csv", d["prism_vols"][:, None])
    save("data_disk_hstar.csv", np.array([[float(d["h_star"])]]))
    save("data_disk_Vprism.csv", np.array([[float(d["V_prism"])]]))
    save("data_disk_Vslices.csv", np.array([[float(d["V_slices"])]]))

    print(f"wrote CSV files to {OUT}")
    print("files:", sorted(os.listdir(OUT)))


if __name__ == "__main__":
    main()
