#!/usr/bin/env wolframscript
(* Rebuild every figure in the post in pure Wolfram Language.

   Input data lives in ./data/ and was produced by export_for_wolfram.py
   in the project root.  We just consume those CSVs here.

   Outputs (next to the notebook):
     sofa.png, sofa.gif                 -- Section 4 hook + Section-2 inset
     area_vs_angle.png                  -- Section 7 curve
     sweep_static.png                   -- (legacy, kept for now)
     sweep.gif                          -- Section 7 alpha-sweep animation
     sofa_3d.gif, gerver_3d.png         -- Section 6 3D scenes
     tilt_attempt.png                   -- Section 6 tilt-failure panel
     aisle_<a>.png, sofa_<a>.png        -- Section 7 grid thumbnails
     thumb_unit_square.png, thumb_half_disk.png,
     thumb_hammersley.png, thumb_gerver.png  -- Table 1 thumbnails
*)

hereDir = DirectoryName[
  If[$ScriptCommandLine =!= {}, $InputFileName, NotebookFileName[]]];
dataDir = FileNameJoin[{hereDir, "data"}];
outDir  = hereDir;

readCSV[name_]  := Import[FileNameJoin[{dataDir, name}], "CSV"];
readVec[name_]  := Flatten @ readCSV[name];

(* ===================================================================
   Reusable styling
   =================================================================== *)

sofaFace = RGBColor[0.81, 0.89, 1.0];
sofaEdge = RGBColor[0.12, 0.31, 0.66];
trajCol  = RGBColor[0.82, 0.41, 0.12];
aisleFace = GrayLevel[0.94];
aisleEdge = GrayLevel[0.4];

(* ===================================================================
   Figure 2 (static optimised 2D sofa)
   =================================================================== *)

makeFigure2[] := Module[{boundary, trajectory, x0, x1, y0, y1, pad = 0.5},
  boundary   = readCSV["boundary_90.csv"];
  trajectory = readCSV["trajectory_90.csv"];
  {x0, x1} = MinMax[boundary[[All, 1]]];
  {y0, y1} = MinMax[boundary[[All, 2]]];
  Export[FileNameJoin[{outDir, "sofa.png"}],
    Graphics[{
      {EdgeForm[Directive[sofaEdge, AbsoluteThickness[2]]],
       sofaFace,
       Polygon[boundary]},
      {trajCol, Dashed, AbsoluteThickness[1.2],
       Line[trajectory]},
      {trajCol, PointSize[0.012],
       Point[trajectory[[{1, -1}]]]}
      },
      PlotRange -> {{x0 - pad, x1 + pad}, {y0 - pad, y1 + pad}},
      AspectRatio -> Automatic, Frame -> True, FrameTicksStyle -> 8,
      PlotLabel -> Style[StringForm[
        "optimised sofa, area = ``", NumberForm[2.21968, {7, 5}]], 11],
      ImageSize -> 720, Background -> White],
    ImageResolution -> 160];
  Print["wrote sofa.png"];
];

(* ===================================================================
   Figure 1 (animated corridor sweep)
   =================================================================== *)

(* Recompute the L-corridor polygon at (theta, {cx, cy}) directly in
   Wolfram \[LongDash] cheaper and cleaner than round-tripping it through CSV. *)
corridorVerts[theta_, {cx_, cy_}, big_: 40.0] := Module[
  {c = Cos[theta], s = Sin[theta], ex, ey},
  ex = {c, -s}; ey = {s, c};
  Table[{cx, cy} + (p[[1]] - 1) ex + (p[[2]] - 1) ey,
    {p, {{0., 0.}, {big, 0.}, {big, 1.}, {1., 1.}, {1., big}, {0., big}}}]];

makeFigure1[nFrames_: 60] := Module[
  {boundary, traj, thetas, frames, frameIdx, clipBox},
  boundary = readCSV["boundary_90.csv"];
  traj     = readCSV["trajectory_90.csv"];
  thetas   = readVec["thetas_90.csv"];
  (* sub-sample to nFrames *)
  frameIdx = Round @ Subdivide[1, Length[thetas], nFrames - 1];
  (* visible window *)
  clipBox = Rectangle[{-2.5, -3.0}, {2.5, 2.0}];
  frames = Table[
    With[{i = frameIdx[[k]],
          corr = corridorVerts[thetas[[#]], traj[[#]]] &},
      Graphics[{
        {EdgeForm[Directive[sofaEdge, AbsoluteThickness[2]]],
         sofaFace, Polygon[boundary]},
        {EdgeForm[Directive[Darker[Red, 0.2], AbsoluteThickness[1]]],
         Opacity[0.40, RGBColor[1.0, 0.94, 0.84]],
         (* clip the (huge) corridor polygon to the visible window so the
            scene doesn't zoom out as the arms swing through *)
         RegionIntersection[Polygon[corridorVerts[thetas[[i]], traj[[i]]]],
                             clipBox]},
        Text[Style[StringForm["\[Theta] = ``\[Degree]",
          NumberForm[N[thetas[[i]] 180/Pi], {4, 1}]], 12],
          Scaled[{0.05, 0.96}], {-1, 1}]},
       PlotRange -> {{-2.5, 2.5}, {-3.0, 2.0}}, AspectRatio -> Automatic,
       ImageSize -> 640, Frame -> False, Background -> White]],
    {k, Length[frameIdx]}];
  Export[FileNameJoin[{outDir, "sofa.gif"}], frames,
    "DisplayDurations" -> 0.06, "AnimationRepetitions" -> Infinity];
  Print["wrote sofa.gif (", Length[frames], " frames)"];
];

(* ===================================================================
   Figure 3 (area vs angle curve)
   =================================================================== *)

makeFigure3[] := Module[{summary, alphas, areas, plot},
  summary = Import[FileNameJoin[{dataDir, "sweep_summary.json"}], "RawJSON"];
  alphas = summary["alphas_deg"]; areas = summary["areas"];
  plot = Show[{
    Plot[2.22 (Pi/2)/(Pi - a Degree), {a, 60, 175},
      PlotStyle -> Directive[Dashed, RGBColor[0.8, 0.55, 0.0]],
      PlotLegends -> Placed[{"heuristic C/(\[Pi]\[Minus]\[Alpha])"}, {0.25, 0.85}]],
    ListLinePlot[Transpose[{alphas, areas}],
      PlotStyle -> Directive[sofaEdge, AbsoluteThickness[2]],
      PlotMarkers -> {Automatic, Medium},
      PlotLegends -> Placed[{"optimised sofa"}, {0.25, 0.92}]],
    Graphics[{LightGray, Dashed, Line[{{50, Pi/2. + 2./Pi}, {175, Pi/2. + 2./Pi}}]}]},
    PlotRange -> {{50, 175}, {0, 6}},
    AxesLabel -> {"bend angle \[Alpha] (\[Degree])", "sofa area S(\[Alpha])"},
    Frame -> True, GridLines -> Automatic, ImageSize -> 720,
    PlotLabel ->
      Style["Sofa area vs corridor bend angle", 12, FontFamily -> "Helvetica"]];
  Export[FileNameJoin[{outDir, "area_vs_angle.png"}], plot,
    ImageResolution -> 160];
  Print["wrote area_vs_angle.png"];
];

(* ===================================================================
   Figure 4 (animated alpha sweep)
   =================================================================== *)

makeFigure4[] := Module[{summary, alphas, areas, frames, win, clipBox},
  summary = Import[FileNameJoin[{dataDir, "sweep_summary.json"}], "RawJSON"];
  alphas = summary["alphas_deg"]; areas = summary["areas"];
  (* Common window across all frames so the animation reads as a single
     scene growing rather than each frame being re-zoomed. *)
  win = {{-4.5, 4.5}, {-4.5, 4.5}};
  clipBox = Rectangle[{win[[1, 1]], win[[2, 1]]}, {win[[1, 2]], win[[2, 2]]}];
  frames = Table[
    Module[{aisle, sofa, clip},
      aisle = readCSV["aisle_" <> ToString[a] <> ".csv"];
      sofa  = readCSV["sofa_"  <> ToString[a] <> ".csv"];
      clip = RegionIntersection[Polygon[aisle], clipBox];
      Graphics[{
        {aisleFace, EdgeForm[Directive[aisleEdge, AbsoluteThickness[0.8]]],
         clip},
        {sofaFace, EdgeForm[Directive[sofaEdge, AbsoluteThickness[1.6]]],
         Polygon[sofa]},
        Text[Style[StringForm["\[Alpha] = ``\[Degree]    area = ``",
          a, NumberForm[areas[[Position[alphas, a][[1, 1]]]], {5, 3}]], 11],
          Scaled[{0.04, 0.95}], {-1, 1}]},
        PlotRange -> win, AspectRatio -> 1,
        ImageSize -> 520, Background -> White]],
    {a, alphas}];
  Export[FileNameJoin[{outDir, "sweep.gif"}], frames,
    "DisplayDurations" -> 0.7, "AnimationRepetitions" -> Infinity];
  Print["wrote sweep.gif (", Length[frames], " frames)"];
];

(* ===================================================================
   Figure 7 (rigid-tilt fail panel)
   =================================================================== *)

makeFigure7[] := Module[{boundary, prism, cx, cy, frame, height = 1.0},
  boundary = readCSV["boundary_90.csv"];
  cx = Mean[boundary[[All, 1]]]; cy = Mean[boundary[[All, 2]]];
  boundary = Transpose[{boundary[[All, 1]] - cx, boundary[[All, 2]] - cy}];

  (* Build the side faces of a vertical prism with bottom=boundary at z=0,
     top=boundary at z=height. *)
  prism[tiltDeg_] := Module[
    {ry, b3, t3, sides},
    ry = RotationMatrix[tiltDeg Degree, {0, 1, 0}];
    b3 = (ry . #) & /@ (Append[#, 0] & /@ boundary);
    t3 = (ry . #) & /@ (Append[#, height] & /@ boundary);
    sides = Table[Polygon[{b3[[i]], b3[[i + 1]], t3[[i + 1]], t3[[i]]}],
      {i, 1, Length[b3] - 1}];
    {b3, t3, sides}];

  frame[tiltDeg_] := Module[{b3, t3, sides, zMax, zMin, fits, col, edge},
    {b3, t3, sides} = prism[tiltDeg];
    zMax = Max[Join[b3[[All, 3]], t3[[All, 3]]]];
    zMin = Min[Join[b3[[All, 3]], t3[[All, 3]]]];
    fits = zMax <= 1.0 + 1.*^-9 && zMin >= 0.0 - 1.*^-9;
    col = If[fits, sofaFace, RGBColor[1.0, 0.75, 0.75]];
    edge = If[fits, sofaEdge, Darker @ Red];
    Graphics3D[{
      EdgeForm[Directive[edge, AbsoluteThickness[0.6]]],
      Opacity[0.85, col],
      sides, Polygon[b3], Polygon[t3]},
      PlotRange -> {{-1.5, 1.5}, {-2.5, 2.5}, {-0.4, 1.4}},
      BoxRatios -> {1, 1.4, 0.6}, ViewPoint -> {2.6, -3.0, 1.6},
      Axes -> False, ImageSize -> 280, Background -> White,
      PlotLabel -> Style[
        StringForm["tilt = ``\[Degree]\n``",
          tiltDeg, If[fits, "fits in the tube", "pokes through ceiling/floor"]],
        9]]];
  Export[FileNameJoin[{outDir, "tilt_attempt.png"}],
    GraphicsRow[frame /@ {0, 10, 25}, Spacings -> 0,
      ImageSize -> 880, Background -> White],
    ImageResolution -> 160];
  Print["wrote tilt_attempt.png"];
];

(* ===================================================================
   Figure 6 / gerver_3d (3D animated sofa in the L-tube)
   =================================================================== *)

(* The unit-square L-tube as a list of side-face polygons. *)
ltubeFaces[arm_: 4.5, h_: 1.0] := Module[
  {(* world L corridor in xy, then extruded to z in [0, h] *)
   lpts, b, t},
  lpts = {{-arm, 0}, {1, 0}, {1, -arm}, {0, -arm}, {0, 1}, {-arm, 1}, {-arm, 0}};
  b = Append[#, 0] & /@ lpts;
  t = Append[#, h] & /@ lpts;
  {b, t, Table[Polygon[{b[[i]], b[[i + 1]], t[[i + 1]], t[[i]]}],
    {i, 1, Length[lpts] - 1}]}];

makeGerver3D[] := Module[
  {boundary, b, t3, b3, sides, cx, cy, tubeB, tubeT, tubeSides, h = 1.0},
  boundary = readCSV["boundary_90.csv"];
  cx = Mean[boundary[[All, 1]]]; cy = Mean[boundary[[All, 2]]];
  boundary = Transpose[{boundary[[All, 1]] - cx, boundary[[All, 2]] - cy}];
  b3 = Append[#, 0] & /@ boundary;
  t3 = Append[#, h] & /@ boundary;
  sides = Table[Polygon[{b3[[i]], b3[[i + 1]], t3[[i + 1]], t3[[i]]}],
    {i, 1, Length[boundary] - 1}];
  {tubeB, tubeT, tubeSides} = ltubeFaces[];
  Export[FileNameJoin[{outDir, "gerver_3d.png"}],
    Graphics3D[{
      {Opacity[0.10, GrayLevel[0.8]],
       EdgeForm[Directive[GrayLevel[0.5], AbsoluteThickness[0.5]]],
       tubeSides, Polygon[tubeB], Polygon[tubeT]},
      {Opacity[0.85, sofaFace],
       EdgeForm[Directive[sofaEdge, AbsoluteThickness[0.6]]],
       sides, Polygon[b3], Polygon[t3]}
    }, PlotRange -> {{-3, 1.5}, {-3, 2}, {-0.1, 1.4}},
      BoxRatios -> {1.5, 1.5, 0.5}, ViewPoint -> {2.8, -2.5, 1.6},
      ImageSize -> 720, Background -> White,
      Axes -> True, AxesLabel -> {"x", "y", "z"}],
    ImageResolution -> 160];
  Print["wrote gerver_3d.png"];
];

makeFigure6[nFrames_: 50] := Module[
  {boundary, traj, ctrl, thetas, frames, cx, cy, b3, t3, sides, tubeB, tubeT, tubeSides, h = 1.0},
  boundary = readCSV["boundary_90.csv"];
  ctrl     = readCSV["controls_90.csv"];
  thetas   = Subdivide[0.01, Pi/2 - 0.01, nFrames - 1];
  (* For the world-frame animation, lift the 2D Gerver motion of Section 4
     into the lab frame: rotation theta about world z, translated so that
     the sofa-frame corner C(theta) sits at world (1, 1, 0.5). *)
  trajFn[ctrlPts_, ths_] := Module[{m = Length[ctrlPts], cth, fx, fy},
    cth = Subdivide[0., Pi/2., m - 1];
    fx = Interpolation[Transpose[{cth, ctrlPts[[All, 1]]}], InterpolationOrder -> 3];
    fy = Interpolation[Transpose[{cth, ctrlPts[[All, 2]]}], InterpolationOrder -> 3];
    Transpose[{fx /@ ths, fy /@ ths}]];

  (* Body in sofa frame: extruded boundary at z in [0, 1]. *)
  b3 = Append[#, 0.] & /@ boundary;
  t3 = Append[#, h] & /@ boundary;
  {tubeB, tubeT, tubeSides} = ltubeFaces[5.5, h];

  frames = Table[
    Module[{theta = thetas[[k]], c2, c3, Rmat, transB, transT, transSides,
            txform},
      c2 = First @ trajFn[ctrl, {theta}];
      c3 = {c2[[1]], c2[[2]], 0.5};
      Rmat = RotationMatrix[theta, {0, 0, 1}];
      txform[p_] := (Rmat . p) + ({1., 1., 0.5} - Rmat . c3);
      transB = txform /@ b3;
      transT = txform /@ t3;
      transSides = Table[Polygon[{txform[b3[[i]]], txform[b3[[i + 1]]],
                                   txform[t3[[i + 1]]], txform[t3[[i]]]}],
        {i, 1, Length[boundary] - 1}];
      Graphics3D[{
        {Opacity[0.10, GrayLevel[0.8]],
         EdgeForm[Directive[GrayLevel[0.5], AbsoluteThickness[0.5]]],
         tubeSides, Polygon[tubeB], Polygon[tubeT]},
        {Opacity[0.85, sofaFace],
         EdgeForm[Directive[sofaEdge, AbsoluteThickness[0.5]]],
         transSides, Polygon[transB], Polygon[transT]}
      },
        PlotRange -> {{-3, 5.5}, {-3, 5.5}, {-0.1, 1.4}},
        BoxRatios -> {1, 1, 0.25}, ViewPoint -> {3.0, -3.0, 1.8},
        ImageSize -> 640, Background -> White, Axes -> False,
        PlotLabel -> Style[StringForm["\[Theta] = ``\[Degree]",
          NumberForm[N[theta 180/Pi], {4, 1}]], 10]]],
    {k, nFrames}];

  Export[FileNameJoin[{outDir, "sofa_3d.gif"}], frames,
    "DisplayDurations" -> 0.07, "AnimationRepetitions" -> Infinity];
  Print["wrote sofa_3d.gif (", Length[frames], " frames)"];
];

(* ===================================================================
   Per-angle thumbnails (aisle_<a>.png, sofa_<a>.png)
   =================================================================== *)

aisleSofaThumbs[] := Module[{alphas, summary},
  summary = Import[FileNameJoin[{dataDir, "sweep_summary.json"}], "RawJSON"];
  alphas = summary["alphas_deg"];
  (* Compute common bounding box for the sofa thumbnails so they're visually comparable. *)
  Do[
    Module[{aisle, sofa, win = 3.5, scale = 4.5, cx, cy, xs, ys},
      aisle = readCSV["aisle_" <> ToString[a] <> ".csv"];
      sofa  = readCSV["sofa_"  <> ToString[a] <> ".csv"];
      (* aisle thumb: window around inner corner (1, 1) *)
      Export[FileNameJoin[{outDir, "aisle_" <> ToString[a] <> ".png"}],
        Graphics[{
          aisleFace, EdgeForm[Directive[aisleEdge, AbsoluteThickness[1.4]]],
          Polygon[aisle],
          {RGBColor[0.8, 0.2, 0.1], PointSize[0.035], Point[{1, 1}]}},
          PlotRange -> {{1 - win, 1 + win}, {1 - win, 1 + win}},
          AspectRatio -> Automatic, ImageSize -> 240,
          Background -> White],
        ImageResolution -> 140];
      (* sofa thumb: centred on its centroid, common scale *)
      cx = Mean[sofa[[All, 1]]]; cy = Mean[sofa[[All, 2]]];
      xs = sofa[[All, 1]] - cx; ys = sofa[[All, 2]] - cy;
      Export[FileNameJoin[{outDir, "sofa_" <> ToString[a] <> ".png"}],
        Graphics[{
          sofaFace, EdgeForm[Directive[sofaEdge, AbsoluteThickness[1.4]]],
          Polygon[Transpose[{xs, ys}]]},
          PlotRange -> {{-scale, scale}, {-scale, scale}},
          AspectRatio -> Automatic, ImageSize -> 240,
          Background -> White],
        ImageResolution -> 140];
    ],
    {a, alphas}];
  Print["wrote 7 aisle + 7 sofa thumbnails"];
];

(* ===================================================================
   Table-1 shape thumbnails
   =================================================================== *)

table1Thumbs[] := Module[{},
  (* unit square *)
  Export[FileNameJoin[{outDir, "thumb_unit_square.png"}],
    Graphics[{sofaFace, EdgeForm[sofaEdge],
       Polygon[{{0, 0}, {1, 0}, {1, 1}, {0, 1}}]},
      PlotRange -> {{-0.3, 1.3}, {-0.3, 1.3}}, AspectRatio -> 1,
      ImageSize -> 140, Background -> White],
    ImageResolution -> 140];
  (* unit half-disk *)
  Module[{th = Subdivide[0., Pi, 60]},
    Export[FileNameJoin[{outDir, "thumb_half_disk.png"}],
      Graphics[{sofaFace, EdgeForm[sofaEdge],
        Polygon[Join[{{1, 0}}, Transpose[{Cos[th], Sin[th]}], {{-1, 0}}]]},
        PlotRange -> {{-1.3, 1.3}, {-0.3, 1.3}}, AspectRatio -> 1,
        ImageSize -> 140, Background -> White],
      ImageResolution -> 140]];
  (* Hammersley sofa: rectangle width 4/pi + two unit quarter-circles + notch *)
  Module[{r = 2./Pi, thL = Subdivide[Pi, Pi/2, 40],
          thR = Subdivide[Pi/2, 0., 40], thN = Subdivide[0., Pi, 40], boundary},
    boundary = Join[
      Transpose[{-r + Cos[thL], Sin[thL]}],
      Transpose[{ r + Cos[thR], Sin[thR]}],
      {{r + 1, 0}},
      Reverse @ Transpose[{r Cos[thN], r Sin[thN]}],
      {{-r - 1, 0}}];
    Export[FileNameJoin[{outDir, "thumb_hammersley.png"}],
      Graphics[{sofaFace, EdgeForm[sofaEdge], Polygon[boundary]},
        PlotRange -> {{-2, 2}, {-0.3, 1.3}}, AspectRatio -> Automatic,
        ImageSize -> 240, Background -> White],
      ImageResolution -> 140]];
  (* Gerver thumbnail: same as sofa_90 *)
  Module[{b = readCSV["boundary_90.csv"], cx, cy},
    cx = Mean[b[[All, 1]]]; cy = Mean[b[[All, 2]]];
    Export[FileNameJoin[{outDir, "thumb_gerver.png"}],
      Graphics[{sofaFace, EdgeForm[sofaEdge],
        Polygon[Transpose[{b[[All, 1]] - cx, b[[All, 2]] - cy}]]},
        PlotRange -> Automatic, AspectRatio -> Automatic,
        ImageSize -> 180, Background -> White],
      ImageResolution -> 140]];
  Print["wrote 4 table-1 thumbs"];
];

(* ===================================================================
   Driver
   =================================================================== *)

Print["building figures from data in ", dataDir];
t0 = AbsoluteTime[];
makeFigure2[];
makeFigure1[];
makeFigure3[];
makeFigure4[];
makeFigure7[];
makeGerver3D[];
makeFigure6[];
aisleSofaThumbs[];
table1Thumbs[];
Print["all figures rebuilt in ", Round[AbsoluteTime[] - t0], " s"];
