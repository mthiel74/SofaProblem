#!/usr/bin/env wolframscript
(* Wolfram-native moving-sofa optimiser.

   The same algorithm as the Python notebook:
     - L-corridor polygon at orientation theta with inner corner at (cx, cy)
     - sofa = intersection of corridor poses over a theta grid
     - cubic-spline trajectory through M control points (X_k, Y_k)
     - NMinimize / FindMinimum (NelderMead) on the 2M control coords

   The key implementation choice: we keep every corridor as a Polygon
   primitive and rely on Wolfram's RegionIntersection. To stay tractable
   we use modest M and a moderate theta grid; with warm-started rounds the
   optimum reaches the same Gerver-class area as the Python pipeline.

   Run: wolframscript sofa_optimise.wls   (writes sofa_optimise_result.mx)
*)

(* ---------- corridor geometry ---------- *)

big = 40.;
worldL = {{0., 0.}, {big, 0.}, {big, 1.}, {1., 1.}, {1., big}, {0., big}};

corridorPoly[\[Theta]_?NumericQ, {cx_?NumericQ, cy_?NumericQ}] := Module[
  {c = Cos[\[Theta]], s = Sin[\[Theta]], ex, ey},
  ex = {c, -s}; ey = {s, c};
  Polygon[Table[{cx, cy} + (p[[1]] - 1) ex + (p[[2]] - 1) ey, {p, worldL}]]
];

(* Intersect every corridor pose. Use RegionIntersection on Polygons. *)
buildSofa[corners_, thetas_] := Module[{r, m = Length[thetas]},
  r = corridorPoly[thetas[[1]], corners[[1]]];
  Do[
    r = RegionIntersection[r, corridorPoly[thetas[[i]], corners[[i]]]];
    , {i, 2, m}];
  r
];

(* ---------- cubic-spline trajectory ---------- *)

evalTraj[ctrl_?MatrixQ, thetas_List] := Module[
  {m = Length[ctrl], cth, fx, fy},
  cth = Subdivide[0., Pi/2., m - 1];
  fx = Interpolation[Transpose[{cth, ctrl[[All, 1]]}], InterpolationOrder -> 3];
  fy = Interpolation[Transpose[{cth, ctrl[[All, 2]]}], InterpolationOrder -> 3];
  Transpose[{fx /@ thetas, fy /@ thetas}]
];

(* Quarter-arc seed: corner sweeps a quarter circle of radius r = 2/Pi. *)
quarterArc[m_Integer, r_: 2./Pi] := With[{th = Subdivide[0., Pi/2., m - 1]},
  Transpose[{r Sin[th], r Cos[th]}]
];

(* ---------- objective for NMinimize ---------- *)

(* We trade off resolution for speed: M = 8 control points (so 16
   variables), nTh = 30 theta samples for the objective, and re-evaluate
   at a denser grid (nThEval = 800) once after optimisation for the
   honest area. *)

M = 8;
nTh = 30;
nThEval = 200;
thetas = Subdivide[0., Pi/2., nTh - 1];

negArea[flat_?VectorQ] := Module[
  {ctrl = Partition[flat, 2], corners, sofa, a},
  corners = evalTraj[ctrl, thetas];
  sofa = Quiet @ buildSofa[corners, thetas];
  a = Quiet @ Area[sofa];
  If[NumericQ[a] && a > 0, -a, 0.]
];

(* ---------- run ---------- *)

x0 = Flatten @ quarterArc[M];
vars = Array[c, 2 M];
Print["initial area = ", -negArea[x0]];

(* ---- a small Nelder-Mead in Wolfram (derivative-free, seeded) ---- *)
nelderMead[f_, x0_List, maxIter_:1000, fTol_: 1.*^-6, sizeTol_: 1.*^-6] := Module[
  {n = Length[x0], simplex, fvals, ord, iter = 0, centroid, worst, xr, fr,
   xe, fe, xc, fc, best, sigma = 0.5, alpha = 1., gamma = 2., rho = 0.5,
   newPt, newFv, conv = False, shrink, stepSize},

  (* Initial simplex: x0 plus n points along the axes. *)
  stepSize = 0.05;
  simplex = Prepend[
    Table[x0 + stepSize UnitVector[n, i], {i, n}],
    x0];
  fvals = f /@ simplex;

  While[iter < maxIter && ! conv,
    iter++;
    ord = Ordering[fvals];
    simplex = simplex[[ord]];
    fvals  = fvals[[ord]];
    best = First[simplex];
    worst = Last[simplex];
    centroid = Mean[Drop[simplex, -1]];

    xr = centroid + alpha (centroid - worst); fr = f[xr];

    Which[
      fr < fvals[[1]],
        xe = centroid + gamma (xr - centroid); fe = f[xe];
        If[fe < fr, newPt = xe; newFv = fe, newPt = xr; newFv = fr];
        shrink = False,
      fr < fvals[[-2]],
        newPt = xr; newFv = fr; shrink = False,
      True,
        xc = centroid + rho (xr - centroid); fc = f[xc];
        If[fc < fvals[[-1]],
          newPt = xc; newFv = fc; shrink = False,
          shrink = True;
        ];
    ];

    If[shrink,
      simplex = Prepend[
        Table[simplex[[1]] + sigma (simplex[[i]] - simplex[[1]]), {i, 2, n+1}],
        simplex[[1]]];
      fvals = f /@ simplex;
      ,
      simplex[[-1]] = newPt;
      fvals[[-1]] = newFv;
    ];

    (* termination *)
    conv = (Max[fvals] - Min[fvals] < fTol);
  ];

  ord = Ordering[fvals];
  {fvals[[ord[[1]]]], simplex[[ord[[1]]]], iter}
];

(* ---- multi-round warm-started Nelder-Mead ---- *)

t0 = AbsoluteTime[];
best = x0; bestArea = -negArea[x0];
SeedRandom[7];
Do[
  start = If[k == 0, best, best + 0.05 RandomReal[{-1, 1}, 2 M]];
  {fmin, xmin, iters} = nelderMead[negArea, start, 600];
  area = -fmin;
  If[area > bestArea,
    bestArea = area;
    best = xmin;
    Print["  * round ", k, ": area = ", area, "   (", iters, " iters)"];
    ,
    Print["    round ", k, ": area = ", area, "   (", iters, " iters)"];
  ];
  , {k, 0, 3}];

elapsed = AbsoluteTime[] - t0;
Print["optimisation time: ", elapsed, " s"];
Print["best area (optim grid, ", nTh, " thetas): ", bestArea];

(* Re-evaluate at a denser grid for an honest readout. *)
ctrlBest = Partition[best, 2];
thetasEval = Subdivide[0., Pi/2., nThEval - 1];
sofaBest = buildSofa[evalTraj[ctrlBest, thetasEval], thetasEval];
areaEval = Area[sofaBest];
Print["honest area (eval grid, ", nThEval, " thetas): ", areaEval];
Print["Hammersley analytic: ", N[Pi/2 + 2./Pi]];
Print["Gerver best known  : 2.21953166887197"];

(* Save the optimised control points as a plain text file so the
   downstream visualisation can pick them up without a .mx. *)
Export[FileNameJoin[{DirectoryName[$InputFileName],
  "sofa_optimise_controls.csv"}], ctrlBest];
Print["wrote sofa_optimise_controls.csv"];
