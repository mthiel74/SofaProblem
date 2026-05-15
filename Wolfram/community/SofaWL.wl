(* ::Package:: *)

(* ::Title:: *)
(*SofaWL: a Wolfram-Language port of the Python moving-sofa optimiser*)

(* ::Text:: *)
(*Literal translation of sofa_pipeline.py (sofa-fixed-frame intersection-of-corridors formulation). The trajectory of the L's inner corner is a cubic-spline interpolation through M control points; the optimiser maximises the area of the intersection polygon over the (X,Y) values at those control points.*)

BeginPackage["SofaWL`"];

CorridorPolygon::usage = "CorridorPolygon[theta, {cx, cy}] returns the L-corridor polygon at orientation theta with inner corner at {cx, cy} in the sofa-fixed frame.";

BuildSofa::usage = "BuildSofa[corners, thetas] returns the polygon intersection of every corridor pose. corners has length Length[thetas].";

EvaluateTrajectory::usage = "EvaluateTrajectory[controlXY, thetas] cubic-spline interpolates M control points (X,Y) at thetas in [0, Pi/2].";

NegArea::usage = "NegArea[flat, thetas, M] is -Area[sofa] for use with NMinimize / FindMinimum. flat is a length-(2 M) real vector.";

QuarterArcControls::usage = "QuarterArcControls[M, r] returns M (X,Y) control points along a quarter circle of radius r.";

OptimiseSofa::usage = "OptimiseSofa[M, nThetas, nRounds, maxIter] runs Nelder-Mead warm-started rounds and returns <|\"Controls\" -> ctrlXY, \"Area\" -> a, \"Sofa\" -> polygon|>.";

Begin["`Private`"];

(* ---------- corridor polygon ------------------------------------ *)

big = 40.;
(* World-frame L vertices, CCW, inner corner at (1,1). *)
worldL = {{0., 0.}, {big, 0.}, {big, 1.}, {1., 1.}, {1., big}, {0., big}};

CorridorPolygon[theta_?NumericQ, {cx_?NumericQ, cy_?NumericQ}] := Module[
  {c = Cos[theta], s = Sin[theta], ex, ey},
  ex = {c, -s}; ey = {s, c};
  Polygon[Table[{cx, cy} + (p[[1]] - 1) ex + (p[[2]] - 1) ey, {p, worldL}]]
];

(* ---------- intersect every pose -------------------------------- *)

BuildSofa[corners_List, thetas_List] := Module[{poly, nxt, i, n = Length[thetas]},
  poly = DiscretizeGraphics @ Graphics @ CorridorPolygon[thetas[[1]], corners[[1]]];
  Do[
    nxt = DiscretizeGraphics @ Graphics @ CorridorPolygon[thetas[[i]], corners[[i]]];
    poly = RegionIntersection[poly, nxt];
    If[Quiet@Area[poly] == 0., Return[poly, Module]],
    {i, 2, n}];
  poly
];

(* ---------- trajectory ----------------------------------------- *)

EvaluateTrajectory[controlXY_?MatrixQ, thetas_List] := Module[
  {m = Length[controlXY], cthetas, fx, fy},
  cthetas = Subdivide[0., Pi/2., m - 1];
  fx = Interpolation[Transpose[{cthetas, controlXY[[All, 1]]}], InterpolationOrder -> 3];
  fy = Interpolation[Transpose[{cthetas, controlXY[[All, 2]]}], InterpolationOrder -> 3];
  Transpose[{fx /@ thetas, fy /@ thetas}]
];

QuarterArcControls[m_Integer, r_: 2./Pi] := Module[{th = Subdivide[0., Pi/2., m - 1]},
  Transpose[{r Sin[th], r Cos[th]}]
];

NegArea[flat_?VectorQ, thetas_List, m_Integer] := Module[
  {ctrl = Partition[flat, 2], corners, sofa},
  corners = EvaluateTrajectory[ctrl, thetas];
  sofa = BuildSofa[corners, thetas];
  If[Head[sofa] === MeshRegion || Head[sofa] === BoundaryMeshRegion, -Quiet@Area[sofa], 0.]
];

(* ---------- optimiser ------------------------------------------ *)

OptimiseSofa[m_: 10, nThetas_: 100, nRounds_: 3, maxIter_: 1500] := Module[
  {thetas, x0, vars, best, area, ctrlXY, sofa, k, res, rng = RandomReal[{-1, 1}, #] &},
  thetas = Subdivide[0., Pi/2., nThetas - 1];
  x0 = Flatten[QuarterArcControls[m]];
  vars = Table[Unique["c"], {2 m}];
  best = x0;
  (* repeated warm starts: NMinimize with Nelder-Mead, seeded from the previous best *)
  Do[
    With[{x = If[k == 0, best, best + 0.05 rng[2 m]]},
      res = FindMinimum[
        NegArea[vars, thetas, m],
        Transpose[{vars, x}],
        Method -> {"NelderMead", "PostProcess" -> False},
        MaxIterations -> maxIter, AccuracyGoal -> 6, PrecisionGoal -> 6];
      If[res[[1]] < NegArea[best, thetas, m],
        best = vars /. res[[2]]
      ]
    ],
    {k, 0, nRounds - 1}];
  ctrlXY = Partition[best, 2];
  sofa = BuildSofa[EvaluateTrajectory[ctrlXY, Subdivide[0., Pi/2., 1599]], Subdivide[0., Pi/2., 1599]];
  area = Quiet@Area[sofa];
  <|"Controls" -> ctrlXY, "Area" -> area, "Sofa" -> sofa, "Thetas" -> thetas|>
];

End[];
EndPackage[];
