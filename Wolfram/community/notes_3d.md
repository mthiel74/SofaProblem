# 3D Extension of the Moving Sofa Problem (research notes for the post)

**Formulation.** The natural 3D analogue replaces the L-shaped planar corridor by an L-shaped *tube*: two semi-infinite prisms with unit cross-section meeting at a right angle. The cross-section is usually taken to be the unit square (matching the 2D unit-width hallway), giving a tube `C = ({(x,y,z): 0 ≤ y,z ≤ 1, x ≤ 1}) ∪ ({0 ≤ x,z ≤ 1, y ≤ 1})`. One asks for the rigid body `S ⊂ R^3` of maximum volume that admits a continuous motion `g(t) ∈ SE(3)` with `g(t) S ⊂ C` for all `t`, starting deep in one arm and ending deep in the other. A disk-shaped cross-section is also natural and is sometimes preferred because it removes the privileged "floor" direction.

**Known results.** The 3D problem is essentially open. Almost all literature treats only the 2D version (Hammersley 1968 with area π/2 + 2/π ≈ 2.2074; Gerver 1992 with A_G ≈ 2.21953; Kallus–Romik 2018 upper bound ≤ 2.37; Baek 2024 proof of Gerver's optimality). The only rigorous 3D result I found is for rectangular boxes through L-shaped tubes (Horváth et al.): the optimal *parallelepiped* in a square-section tube is the trivial extruded 2D rectangle. No non-trivial 3D construction beating A_G × 1 ≈ 2.21953 appears in the published literature, although it is widely *expected* that one exists.

**Trivial bounds.**
- Lower bound: V_3 ≥ A_G ≈ 2.21953 (extrude Gerver's sofa along the z-axis).
- Upper bound: no published explicit value; folklore ≲ 3 from swept-volume arguments.

**Why 3D is genuinely harder (not just bigger).** Above the 2D motion you gain two additional rotational axes (pitch and roll about the direction of travel). The sofa can tilt and corkscrew as it rounds the corner, lifting its "tail" into the ceiling clearance of the *other* arm. This breaks the up–down symmetry of the extruded solution and almost certainly permits a strictly larger body. Configuration space jumps from SE(2) (3-DOF) to SE(3) (6-DOF), and rotations no longer commute, so Gerver's local-optimality ODE machinery does not transfer cleanly.

**Algorithmic generalisation.** Your "intersect rotated/translated corridors in the sofa-fixed frame" formulation does extend: discretise a path g(t) = (R(t), p(t)) with R(t) ∈ SO(3) parameterised by, e.g., a quaternion spline or three Euler-angle splines, then S = ∩_t g(t)^-1 C. Optimise control points to maximise vol(S). The "inner edge" of the corner is now a line, not a point, so the natural constraint is that this edge slides along a curve on ∂S. Expect ∼50–200 parameters and millions of voxel/half-space evaluations per objective call; a serious workstation run is feasible (hours to days), but tuning to avoid local minima makes this research-grade rather than turn-the-crank.

### Bounds table

| Setting | Lower bound | Upper bound | Source |
|---|---|---|---|
| 2D, 90° (Gerver sofa) | 2.21953 | 2.21953 (proved optimal) | Gerver 1992; Baek 2024 |
| 2D, 90° earlier | 2.2074 (Hammersley) | 2.37 | Hammersley 1968; Kallus–Romik 2018 |
| 3D, 90° square tube | ≥ 2.21953 (extrusion) | no published explicit value; ≲ 3 folklore | this work / folklore |
| 3D, 90° circular tube | ≥ π/4 × (some 2D-disk-sofa) | open | open |

### Sources
- J. L. Gerver, "On moving a sofa around a corner," *Geom. Dedicata* 42 (1992).
- Y. Kallus, D. Romik, "Improved upper bounds in the moving sofa problem," *Adv. Math.* 340 (2018), arXiv:1706.06630.
- J. Baek, "Optimality of Gerver's Sofa," arXiv:2411.19826 (2024).
- H. T. Croft, K. J. Falconer, R. K. Guy, *Unsolved Problems in Geometry*, Springer 1991 (problem G5).
- Á. G. Horváth et al., "Moving rectangular sofas in planar and spatial corridors."
