#!/usr/bin/env wolframscript
(* Build all figures for the second notebook (5 explorations).
   Reads Python-produced .npz files via numpy text dumps that we make
   along the way. *)

hereDir = DirectoryName[
  If[$ScriptCommandLine =!= {}, $InputFileName, NotebookFileName[]]];

sofaFace = RGBColor[0.81, 0.89, 1.0];
sofaEdge = RGBColor[0.12, 0.31, 0.66];
trajCol  = RGBColor[0.82, 0.41, 0.12];
aisleFace = GrayLevel[0.94];
aisleEdge = GrayLevel[0.4];
warnCol  = RGBColor[0.85, 0.25, 0.20];

readCSV[name_] := Import[FileNameJoin[{hereDir, name}], "CSV"];

(* Common Graphics options *)
basePlotOpts = Sequence[Frame -> True, GridLines -> Automatic,
                         ImageSize -> 640, Background -> White,
                         LabelStyle -> {12, FontFamily -> "Helvetica"}];

(* ===================================================================
   FIG #4: variable-width A(w) curve
   =================================================================== *)

makeWidthCurve[] := Module[{widths, areas, plot, ws, lin},
  widths = readCSV["data_widths.csv"][[All, 1]];
  areas  = readCSV["data_areas.csv"][[All, 1]];
  ws = N[Range[0.0, 1.0, 0.01]];
  (* Simple linear fit A(w) ~ slope * w + intercept *)
  lin = Fit[Transpose[{widths, areas}], {1, w}, w];
  plot = Show[{
    Plot[lin, {w, 0, 1}, PlotStyle -> {Dashed, Gray},
         PlotLegends -> Placed[{"linear fit"}, {0.7, 0.15}]],
    ListPlot[Transpose[{widths, areas}],
             PlotStyle -> Directive[sofaEdge, AbsoluteThickness[2], PointSize[0.012]],
             Joined -> True, PlotLegends -> Placed[{"optimised"}, {0.7, 0.25}]],
    Graphics[{Dashed, Gray,
              Line[{{0, Pi/2 + 2/Pi}, {1, Pi/2 + 2/Pi}}],
              Text[Style["Gerver = 2.2195", 10], {0.55, 2.30}]}]},
    basePlotOpts,
    FrameLabel -> {"arm-2 width  w  (arm-1 width = 1)",
                   "max 2D sofa area  A(w)"},
    PlotRange -> {{0, 1.05}, {0, 2.5}},
    PlotLabel -> Style["Sofa area in an asymmetric L-corridor", 13]];
  Export[FileNameJoin[{hereDir, "fig_width_curve.png"}], plot,
         ImageResolution -> 160];
  Print["wrote fig_width_curve.png"];
];

(* ===================================================================
   FIG #2: disk-tube volume vs slab height
   =================================================================== *)

makeDiskTubePlot[] := Module[{hs, vs, hstar, Vprism, Vslice, plot},
  hs     = readCSV["data_disk_hs.csv"][[All, 1]];
  vs     = readCSV["data_disk_vs.csv"][[All, 1]];
  hstar  = readCSV["data_disk_hstar.csv"][[1, 1]];
  Vprism = readCSV["data_disk_Vprism.csv"][[1, 1]];
  Vslice = readCSV["data_disk_Vslices.csv"][[1, 1]];
  plot = Show[{
    ListLinePlot[Transpose[{hs, vs}],
       PlotStyle -> Directive[sofaEdge, AbsoluteThickness[2]],
       PlotLegends -> Placed[{"V_prism(h) = h \[CenterDot] A(\[Sqrt](1-h\.b2))"}, {0.65, 0.85}]],
    Graphics[{Dashed, Gray,
              Line[{{0, Vslice}, {1, Vslice}}],
              Text[Style[Row[{"slice-upper-bound \[Integral]A(w(z))dz = ",
                              NumberForm[Vslice, {5, 4}]}], 10], {0.6, Vslice + 0.06}],
              Line[{{0, 2.21953}, {1, 2.21953}}],
              Text[Style["square-tube extrusion = 2.2195", 10], {0.5, 2.25}],
              Red, PointSize[0.014], Point[{hstar, Vprism}],
              Text[Style[Row[{"h\.b9 = ", NumberForm[hstar, {3, 3}],
                              ", V = ", NumberForm[Vprism, {4, 3}]}], 10],
                   {hstar, Vprism - 0.18}]}]},
    basePlotOpts,
    FrameLabel -> {"prism height h", "sofa volume V"},
    PlotRange -> {{0, 1.05}, {0, 2.5}},
    PlotLabel -> Style["Disk-cross-section 3D L-tube: prismatic vs slice bound", 13]];
  Export[FileNameJoin[{hereDir, "fig_disk_tube.png"}], plot,
         ImageResolution -> 160];
  Print["wrote fig_disk_tube.png"];
];

(* ===================================================================
   FIG #3: S-corridor sofa
   =================================================================== *)

makeSCorridor[] := Module[{boundary, plot},
  boundary = readCSV["data_s_boundary.csv"];
  plot = Graphics[{
    {EdgeForm[Directive[sofaEdge, AbsoluteThickness[2]]],
     sofaFace, Polygon[boundary]}},
    AspectRatio -> Automatic, ImageSize -> 720,
    Frame -> True, FrameLabel -> {"x", "y"},
    PlotLabel -> Style[
      "S-corridor sofa (two consecutive 90\[Degree] bends)  \[LongDash]  area = 1.013", 13]];
  Export[FileNameJoin[{hereDir, "fig_s_corridor.png"}], plot,
         ImageResolution -> 160];
  Print["wrote fig_s_corridor.png"];
];

(* ===================================================================
   FIG #5: Convex sofa shape vs Gerver shape
   =================================================================== *)

makeConvex[] := Module[{convexBoundary, gerverBoundary, areaC},
  convexBoundary = readCSV["data_convex_boundary.csv"];
  gerverBoundary = readCSV["data_gerver_boundary.csv"];
  (* centre both for side-by-side *)
  Module[{cx, cy},
    cx = Mean[gerverBoundary[[All, 1]]]; cy = Mean[gerverBoundary[[All, 2]]];
    gerverBoundary = Transpose[{gerverBoundary[[All, 1]] - cx,
                                gerverBoundary[[All, 2]] - cy}]];
  Module[{cx, cy},
    cx = Mean[convexBoundary[[All, 1]]]; cy = Mean[convexBoundary[[All, 2]]];
    convexBoundary = Transpose[{convexBoundary[[All, 1]] - cx,
                                 convexBoundary[[All, 2]] - cy}]];
  Export[FileNameJoin[{hereDir, "fig_convex.png"}],
    GraphicsRow[{
      Graphics[{EdgeForm[Directive[sofaEdge, AbsoluteThickness[2]]],
                sofaFace, Polygon[gerverBoundary]},
               PlotRange -> {{-1.5, 1.5}, {-2, 2}}, AspectRatio -> Automatic,
               PlotLabel -> Style["Gerver-class\narea \[TildeTilde] 2.22", 12]],
      Graphics[{EdgeForm[Directive[warnCol, AbsoluteThickness[2]]],
                RGBColor[1.0, 0.92, 0.85], Polygon[convexBoundary]},
               PlotRange -> {{-1.5, 1.5}, {-2, 2}}, AspectRatio -> Automatic,
               PlotLabel -> Style["best convex sofa found\narea \[TildeTilde] 0.81", 12]]},
      Spacings -> 0.6, ImageSize -> 720, Background -> White],
    ImageResolution -> 160];
  Print["wrote fig_convex.png"];
];

(* ===================================================================
   FIG #4 (b): four sofa shapes at varying widths, side by side
   =================================================================== *)

makeWidthGallery[] := Module[{widths, boundaries, panels, p},
  widths = readCSV["data_widths.csv"][[All, 1]];
  panels = Table[
    Module[{b, w = widths[[i]], cx, cy},
      b = readCSV["data_width_" <> StringReplace[ToString[N[w]], "." -> "p"] <> ".csv"];
      cx = Mean[b[[All, 1]]]; cy = Mean[b[[All, 2]]];
      b = Transpose[{b[[All, 1]] - cx, b[[All, 2]] - cy}];
      Graphics[{EdgeForm[Directive[sofaEdge, AbsoluteThickness[1.5]]],
                sofaFace, Polygon[b]},
               PlotRange -> {{-1.5, 1.5}, {-2, 2}}, AspectRatio -> Automatic,
               PlotLabel -> Style[
                 StringForm["w = ``\narea = ``", w,
                            NumberForm[N[readCSV["data_areas.csv"][[i, 1]]], {4, 3}]],
                 10]]],
    {i, Length[widths]}];
  Export[FileNameJoin[{hereDir, "fig_width_gallery.png"}],
    GraphicsRow[panels, Spacings -> 0.3, ImageSize -> 1100, Background -> White],
    ImageResolution -> 130];
  Print["wrote fig_width_gallery.png"];
];

(* ===================================================================
   driver
   =================================================================== *)

t0 = AbsoluteTime[];
makeWidthCurve[];
makeDiskTubePlot[];
makeSCorridor[];
makeConvex[];
makeWidthGallery[];
Print["figures built in ", Round[AbsoluteTime[] - t0], " s"];
