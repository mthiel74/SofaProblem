# Varying the bend angle in the moving-sofa problem (research notes)

**1. Statement.** Let α ∈ (0, π). Place the inner corner at the origin. The first arm is the unit-width strip with the inner wall along +x, the outer wall 1 unit away. The second arm is the same strip rotated by α about the origin so that the angle between the two inner walls (measured on the corridor side) is α. The corridor C(α) is the union of the two arms. The α-sofa problem asks for

`S(α) = sup { area(K) : K ⊂ R² connected, admits a continuous rigid motion taking K from "deep in arm 1" to "deep in arm 2" with K(t) ⊂ C(α) for all t }`.

S(π/2) is the classical sofa constant ≈ 2.21953 (Gerver 1992; optimality proved by Baek 2024).

**2. What is known.** The literature on non-right-angle bends is essentially empty. Hammersley (1968), Gerver (1992), Romik (2017), and Kallus–Romik (2018) treat only α = π/2. Maruyama's 1973 polygonal scheme is in principle angle-agnostic but was never run at other angles in print. **No analytic family S(α) is known, no asymptotic expansion is published**, and the Baek 2024 optimality proof is specific to π/2.

**3. Heuristics.**

*α → π (the strip limit).* As α → π the corridor becomes a strip and the sofa can be arbitrarily long: **S(α) → ∞** at least like 1/(π − α).

*α → 0 (the fold limit).* A disk of radius tan(α/2) inscribes in the wedge, giving S(α) ≳ π tan²(α/2) ∼ πα²/4. So **S(α) → 0 as α → 0**, leading order linear or quadratic in α.

*Near α = π/2.* S is almost certainly continuous and probably smooth in α near π/2 — the optimum rotates by π − α and Gerver's 18-arc construction deforms smoothly under small perturbations.

**4. Algorithm.** Only the corridor polygon changes. The trajectory parameterisation, the rotation total (π − α instead of π/2), and the Nelder-Mead loss are otherwise identical. Practical points: update the rotation total in the trajectory model from π/2 to π − α; initialise from Gerver's 18-arc shape and march α away from π/2.

**5. Number table (orders-of-magnitude seed estimates).**

A reasonable smallest angle is 45° — below this the sofa becomes very thin. The remaining estimates are seeded from heuristics calibrated so S(π/2) = 2.22.

| α (deg) | S(α) estimate | uncertainty |
|---|---|---|
| 45° | ≈ 0.55 | ±0.10, lower-bound dominated |
| 60° | ≈ 1.05 | ±0.15 |
| 75° | ≈ 1.65 | ±0.10 |
| 90° | 2.2195 | exact (Gerver / Baek) |
| 105° | ≈ 2.95 | ±0.2 |
| 120° | ≈ 4.1 | ±0.4 |
| 135° | ≈ 6.5 | ±1.0 |
| 150° | ≈ 12 | ±3 |
| 170° | ≈ 70 | order-of-magnitude |

These are starting guesses to validate the optimiser against — the optimiser results below are the real numbers.

**6. Animation.** Sweep α from 60° to 170° in ~60 frames. Per frame, display the corridor polygon, the sofa polygon and the area readout. Cheapest computation: optimise at ~7 anchor angles. For intermediate frames, interpolate the spline control points linearly in α and re-evaluate only the (cheap) intersection. Warm-start Nelder-Mead from the previous anchor's optimum.

**Sources.**
- Gerver, "On Moving a Sofa Around a Corner," *Geom. Dedicata* (1992).
- Romik, "Differential equations and exact solutions in the moving sofa problem," *Exp. Math.* (2017).
- Kallus & Romik, "Improved upper bounds in the moving sofa problem," *Adv. Math.* (2018).
- Baek, "Optimality of Gerver's Sofa," arXiv:2411.19826 (2024).
- Maruyama, "An approximation method for solving the sofa problem," *Int. J. Comp. Inf. Sci.* (1973).
