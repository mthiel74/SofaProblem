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

(* The unit-square L-tube as a list of side-face polygons, in the SAME
   convention as the 2D optimiser: inner corner at (1, 1), outer corner
   at (0, 0), arms extending in +x and +y from the corner. *)
ltubeFaces[arm_: 4.5, h_: 1.0] := Module[{lpts, b, t},
  lpts = {{0., 0.}, {arm, 0.}, {arm, 1.}, {1., 1.}, {1., arm}, {0., arm},
          {0., 0.}};
  b = Append[#, 0.] & /@ lpts;
  t = Append[#, h]  & /@ lpts;
  {b, t, Table[Polygon[{b[[i]], b[[i + 1]], t[[i + 1]], t[[i]]}],
               {i, 1, Length[lpts] - 1}]}];

(* Cubic spline trajectory like in build_notebook (kept local). *)
evalCorner[ctrl_?MatrixQ, ths_List] := Module[
  {m = Length[ctrl], cth, fx, fy},
  cth = Subdivide[0., Pi/2., m - 1];
  fx = Interpolation[Transpose[{cth, ctrl[[All, 1]]}], InterpolationOrder -> 3];
  fy = Interpolation[Transpose[{cth, ctrl[[All, 2]]}], InterpolationOrder -> 3];
  Transpose[{fx /@ ths, fy /@ ths}]];

(* Lab-frame rigid motion: world = R p + t with R = Rz(theta) Ry(beta)
   chosen so that the sofa-frame corner C(theta, 0.5) lands on the
   world inner corner (1, 1, 0.5). *)
labFrame[theta_, beta_, ctrl_] := Module[{Rz, Ry, R, c2, c3, t},
  Rz = RotationMatrix[theta, {0, 0, 1}];
  Ry = RotationMatrix[beta,  {0, 1, 0}];
  R  = Rz . Ry;
  c2 = First @ evalCorner[ctrl, {theta}];
  c3 = {c2[[1]], c2[[2]], 0.5};
  t  = {1., 1., 0.5} - R . c3;
  {R, t}];

(* Compute the sofa-frame theta at which the sofa is most centred at
   the bend (i.e. (sofa-frame centroid -> world (1, 1)) is solved by
   the labFrame translation). We just pick the theta whose C(theta) is
   closest to the 2D-sofa centroid; that puts the centroid right at the
   world inner corner and the rendering is naturally centred. *)
bestCenteringTheta[boundary_?MatrixQ, ctrl_?MatrixQ] := Module[
  {centroid, thetas, corners, dists},
  centroid = Mean[boundary];
  thetas = Subdivide[0., Pi/2., 200];
  corners = evalCorner[ctrl, thetas];
  dists = Norm /@ ((# - centroid) & /@ corners);
  thetas[[First @ Ordering[dists, 1]]]
];

makeGerver3D[] := Module[
  {boundary, ctrl, R, tvec, b3, t3, sides, tubeB, tubeT, tubeSides, h = 1.0,
   apply, thetaMid, arm = 3.5, cx, cy},
  boundary = readCSV["boundary_90.csv"];
  ctrl     = readCSV["controls_90.csv"];
  thetaMid = bestCenteringTheta[boundary, ctrl];

  (* Sofa in 3D, sofa-fixed frame: extruded boundary, z in [0, 1]. *)
  b3 = Append[#, 0.] & /@ boundary;
  t3 = Append[#, h]  & /@ boundary;
  (* Lab-frame transform: sofa-frame C(theta) -> world (1, 1, 0.5). *)
  {R, tvec} = labFrame[thetaMid, 0., ctrl];
  apply[p_List] := R . p + tvec;
  b3 = apply /@ b3; t3 = apply /@ t3;
  sides = Table[Polygon[{b3[[i]], b3[[i + 1]], t3[[i + 1]], t3[[i]]}],
                {i, 1, Length[b3] - 1}];

  (* Shorter-arm L-tube so the bend dominates the view. *)
  {tubeB, tubeT, tubeSides} = ltubeFaces[arm, h];

  (* Centre the camera on the sofa's projected centroid. *)
  cx = Mean[b3[[All, 1]]]; cy = Mean[b3[[All, 2]]];
  Export[FileNameJoin[{outDir, "gerver_3d.png"}],
    Graphics3D[{
      {Opacity[0.16, GrayLevel[0.85]],
       EdgeForm[Directive[GrayLevel[0.4], AbsoluteThickness[0.6]]],
       tubeSides, Polygon[tubeB], Polygon[tubeT]},
      {Opacity[0.92, sofaFace],
       EdgeForm[Directive[sofaEdge, AbsoluteThickness[0.6]]],
       sides, Polygon[b3], Polygon[t3]}},
     PlotRange -> {{cx - 2.5, cx + 2.5}, {cy - 2.5, cy + 2.5}, {-0.2, 1.4}},
     BoxRatios -> {1, 1, 0.3}, ViewPoint -> {0., 0., 5.0},
     ViewVertical -> {0, 1, 0},
     ImageSize -> 720, Background -> White,
     Boxed -> False, Axes -> False,
     PlotLabel -> Style[
       "Extruded Gerver sofa centred at the bend (theta \[TildeTilde] pi/4); 2D area = 2.2197", 11]],
    ImageResolution -> 160];
  Print["wrote gerver_3d.png"];
];

makeFigure6[nFrames_: 50] := Module[
  {boundary, ctrl, thetas, frames, b3Base, t3Base, tubeB, tubeT, tubeSides,
   h = 1.0, arm = 4.0, cx, cy},
  boundary = readCSV["boundary_90.csv"];
  ctrl     = readCSV["controls_90.csv"];
  thetas   = Subdivide[0.01, Pi/2 - 0.01, nFrames - 1];
  b3Base = Append[#, 0.] & /@ boundary;
  t3Base = Append[#, h]  & /@ boundary;
  {tubeB, tubeT, tubeSides} = ltubeFaces[arm, h];

  (* Compute a centre point for the whole animation: world inner corner
     (1, 1) is the natural centre. *)
  cx = 1.0; cy = 1.0;

  frames = Table[
    Module[{theta = thetas[[k]], R, tvec, apply, b3, t3, sides},
      {R, tvec} = labFrame[theta, 0., ctrl];
      apply[p_List] := R . p + tvec;
      b3 = apply /@ b3Base; t3 = apply /@ t3Base;
      sides = Table[Polygon[{b3[[i]], b3[[i + 1]], t3[[i + 1]], t3[[i]]}],
        {i, 1, Length[b3] - 1}];
      Graphics3D[{
        {Opacity[0.16, GrayLevel[0.85]],
         EdgeForm[Directive[GrayLevel[0.4], AbsoluteThickness[0.6]]],
         tubeSides, Polygon[tubeB], Polygon[tubeT]},
        {Opacity[0.92, sofaFace],
         EdgeForm[Directive[sofaEdge, AbsoluteThickness[0.5]]],
         sides, Polygon[b3], Polygon[t3]}},
        PlotRange -> {{cx - 3.0, cx + 3.0}, {cy - 3.0, cy + 3.0}, {-0.2, 1.4}},
        BoxRatios -> {1, 1, 0.27}, ViewPoint -> {0., 0., 5.0},
        ViewVertical -> {0, 1, 0},
        ImageSize -> 640, Background -> White,
        Boxed -> False, Axes -> False,
        PlotLabel -> Style[StringForm["\[Theta] = ``\[Degree]",
          NumberForm[N[theta 180/Pi], {4, 1}]], 11]]],
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
