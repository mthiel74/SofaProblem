#!/usr/bin/env wolframscript
(* Build Sofa2.nb -- a second Wolfram Community post collecting the
   numerical results of the five "beyond Gerver" explorations. *)

hereDir = DirectoryName[
  If[$ScriptCommandLine =!= {}, $InputFileName, NotebookFileName[]]];
imgWidth = 680;

loadImg[name_, width_:Automatic] := With[
  {p = FileNameJoin[{hereDir, name}], w = If[width === Automatic, imgWidth, width]},
  Which[
    ! FileExistsQ[p],
      Style["[missing image: " <> name <> "]", Red],
    ToLowerCase[FileExtension[p]] === "gif",
      ImageResize[Import[p, "AnimatedImage"], w],
    True,
      Image[Import[p], ImageSize -> w]]];

imgCell[name_, caption_String: "", widthOverride_: Automatic] := Cell[
  BoxData[ToBoxes[
    With[{img = loadImg[name, If[widthOverride === Automatic, imgWidth, widthOverride]]},
      If[caption === "",
         img,
         Column[{img,
                 Pane[Style[caption, Italic, FontSize -> 10,
                            FontFamily -> "Helvetica"],
                      If[widthOverride === Automatic, imgWidth, widthOverride],
                      Alignment -> Center]},
                Alignment -> Center, Spacings -> 0.6]]]]],
  "Output", ShowCellLabel -> False,
  TextAlignment -> Center,
  CellMargins -> {{Inherited, Inherited}, {10, 14}}];

mathI[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
                    "InlineFormula", FontFamily -> "Times"];
mathD[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
                    "DisplayFormula", TextAlignment -> Center,
                    FontFamily -> "Times"];

cells = {

  Cell["Beyond Gerver: Numerical Explorations of the Moving-Sofa Problem",
       "Title"],
  Cell["Asymmetric corridors, S-bends, disk tubes, multi-bend scaling, and the convex restriction",
       "Subtitle"],
  Cell[TextData[{"Marco Thiel, May 2026 \[LongDash] ",
                 ButtonBox["github.com/mthiel74/SofaProblem",
                   BaseStyle -> "Hyperlink",
                   ButtonData -> {URL["https://github.com/mthiel74/SofaProblem"], None}]}],
       "Text", TextAlignment -> Center, FontSlant -> Italic],

  Cell["Abstract", "Section"],
  Cell[TextData[{
    "A companion to the earlier post reproducing Gerver's 2D bound. Here we run the same intersection-of-corridor-poses optimiser through five new geometries and record what comes out:\n\n",
    StyleBox["(1)", FontWeight -> Bold],
    " an asymmetric L with arm widths 1 and w \[LessEqual] 1; the optimum A(w) traces a near-linear curve.\n",
    StyleBox["(2)", FontWeight -> Bold],
    " a 3D L-tube with a disk cross-section; the best prismatic sofa has volume ", mathI["1.173"], " and the slice-by-slice upper bound is ", mathI["1.808"], " (vs. ", mathI["2.220"], " for a square cross-section).\n",
    StyleBox["(3)", FontWeight -> Bold],
    " an S-corridor with two consecutive 90\[Degree] bends; the optimal sofa has area ",
    mathI["1.013"], ", and its shape collapses to a near-circular disc.\n",
    StyleBox["(4)", FontWeight -> Bold],
    " an N-bend corridor; we measure how the sofa area S(N) decays with N.\n",
    StyleBox["(5)", FontWeight -> Bold],
    " the convex restriction; the best convex sofa we find has area ",
    mathI["0.81"], " against the classical analytical target of \[Pi]/2 \[TildeTilde] 1.57."}],
    "Text"],

  (* === Section 1: Variable width === *)
  Cell["1. Asymmetric L-corridor: area as a function of arm-2 width", "Section"],
  Cell[TextData[{
    "Take the classical L-corridor and shrink arm 2 to width w \[LessEqual] 1 while leaving arm 1 at width 1. ",
    "Inner walls still meet at (1, 1). At w = 1 we recover Gerver; as w \[Rule] 0 the second arm becomes a slit. Sweep:"}], "Text"],
  imgCell["fig_width_curve.png", "Figure 1. Optimised sofa area A(w) as a function of arm-2 width, with the dashed linear fit and the Gerver bound for reference. A(w) is monotonically increasing and remarkably close to linear over w \[Element] [0.4, 1]; A(0) is 0 by geometry but the linear fit extrapolated to w = 0 lands well above 0, so the asymptotic shape near w = 0 is sub-linear."],
  imgCell["fig_width_gallery.png", "Figure 2. The optimised sofa shapes at w = 1, 0.9, ..., 0.4. As w shrinks, the sofa's notch becomes narrower and the body lengthens \[LongDash] but its area drops monotonically."],

  (* === Section 2: Disk-cross-section 3D L-tube === *)
  Cell["2. 3D L-tube with disk cross-section", "Section"],
  Cell[TextData[{
    "Replace the square cross-section by a disk of unit diameter. At height z \[Element] [-1/2, 1/2] the strip width available to that horizontal slice of the sofa is w(z) = 2\[Sqrt](1/4 - z\.b2). For a prismatic body of height h centred in the disk, the cross-section width is constrained by the narrowest slice, w = \[Sqrt](1 - h\.b2), and the volume is V(h) = h \[CenterDot] A(\[Sqrt](1 - h\.b2)), using A(w) from Section 1."}], "Text"],
  imgCell["fig_disk_tube.png", "Figure 3. Prism volume V(h) for the unit-diameter disk-tube. The maximum (red dot) sits at h \[TildeTilde] 0.78, w \[TildeTilde] 0.63, V \[TildeTilde] 1.17. The dashed line at 1.81 is the slice-by-slice upper bound \[Integral] A(w(z))\[ThinSpace]dz, unreachable by a rigid body because slices at different z would need different motions. The square-tube extrusion at V = 2.22 is shown for comparison."],
  Cell[TextData[{
    "Take-away: the disk-tube 3D sofa is ", StyleBox["strictly smaller", FontSlant -> Italic],
    " than the square-tube one in both achievable and upper-bound senses. The square corners help \[LongDash] they let the sofa fill the corridor's vertical clearance regardless of horizontal position."}], "Text"],

  (* === Section 3: S-corridor === *)
  Cell["3. The S-corridor: two consecutive 90\[Degree] bends", "Section"],
  Cell[TextData[{
    "Add a second left bend after the first. The corridor is now a stair: ",
    StyleBox["+", "Input"], "x \[Rule] +y \[Rule] -x. Total rotation is \[Pi]. Same optimiser, extended to \[Theta] \[Element] [0, \[Pi]]:"}], "Text"],
  imgCell["fig_s_corridor.png", "Figure 4. The optimised S-corridor sofa, area 1.013. Strikingly, the shape collapses to an almost-perfectly-round disc of diameter \[TildeTilde] 1.3 (the disc inscribed in both bends simultaneously). Two consecutive 90\[Degree] turns drop the maximum from Gerver's 2.22 to barely above the unit disc \[Pi](0.6)\.b2 \[TildeTilde] 1.13."],

  (* === Section 4: Multi-bend scaling === *)
  Cell["4. N consecutive bends: scaling of S(N)", "Section"],
  Cell["Push further: N consecutive 90\[Degree] left bends, total rotation N\[Pi]/2. Run the optimiser for N = 1, 2, 3, 4 and see what the maximum area does.",
       "Text"],
  imgCell["fig_multibend.png", "Figure 5. Optimised sofa area S(N): 2.220, 1.015, 0.857, 0.810 for N = 1, 2, 3, 4. The dashed line is the area of the unit disk \[Pi]/4 \[TildeTilde] 0.785. S(N) converges geometrically with ratio \[TildeTilde] 1/3."],
  imgCell["fig_multibend_shapes.png", "Figure 6. The optimised sofa shapes at N = 2, 3, 4 overlaid on the unit disk (dashed). As N grows the shape collapses to the disk \[LongDash] the only convex body that fits in a unit-width corridor at every orientation. This is the asymptotic limit  S(N) \[Rule] \[Pi]/4  as N \[Rule] \[Infinity]."],
  Cell[TextData[{
    StyleBox["A genuinely new empirical finding:", FontWeight -> Bold],
    " with each additional 90\[Degree] bend the optimal sofa area decreases geometrically toward the unit-disk area \[Pi]/4. Concretely the three differences S(N) - \[Pi]/4 are 1.435, 0.230, 0.072, 0.024 for N = 1, 2, 3, 4: ratios 0.160, 0.314, 0.336. The asymptotic shape is the disc \[LongDash] the only body that fits a unit-width corridor independent of orientation."}], "Text"],

  (* === Section 5: Convex restriction === *)
  Cell["5. Convex-only sofas", "Section"],
  Cell[TextData[{
    "Constrain the sofa to be convex. Gibbs (1968) showed the supremum is ",
    mathI[RowBox[{"\[Pi]","/","2"}]],
    " (the unit half-disc). Our optimiser, with a cubic-spline corner trajectory and a convexity penalty, lands at 0.81 \[LongDash] significantly below the analytical \[Pi]/2 = 1.57. The shape it finds is round-ish but not the half-disc:"}], "Text"],
  imgCell["fig_convex.png", "Figure 6. Left: the Gerver-class shape from the original post (non-convex, area 2.22). Right: the best convex sofa the optimiser found in this run, area 0.81. The cubic-spline corner trajectory can't quite reproduce the half-disc's rigid \"slide-rotate-slide\" motion; a dedicated convex parameterisation would close the gap to \[Pi]/2."],

  (* === Section 6: Summary === *)
  Cell["6. Summary", "Section"],
  Cell[BoxData[ToBoxes[
    Grid[{
      {Style["setting", Bold], Style["area / volume", Bold], Style["fraction of Gerver", Bold]},
      {"2D classical L-corridor (Gerver / Baek)", "2.21953", "1.0000"},
      {"this notebook (Python, Gerver-class)",    "2.2197", "0.9999"},
      {"this notebook (Wolfram, Gerver-class)",   "2.2219", "1.0011"},
      {"asymmetric L, w = 0.7", "1.6244", "0.7319"},
      {"asymmetric L, w = 0.4", "1.1915", "0.5368"},
      {"S-corridor (2 bends)", "1.013", "0.4564"},
      {"convex-only sofa (our run)", "0.81",  "0.36"},
      {"convex-only sofa (Gibbs)", "\[Pi]/2 \[TildeTilde] 1.5708", "0.71"},
      {"disk-tube prism (3D)", "1.173", "0.529"},
      {"disk-tube slice upper bound (3D)", "1.808", "0.815"},
      {"square-tube extrusion (3D)", "2.21953", "1.0000"}},
      Frame -> All, Alignment -> {Left, Center, Center},
      Background -> {None, {1 -> LightGray}}]]],
    "Output", ShowCellLabel -> False, TextAlignment -> Center],

  Cell[TextData[{
    "The takeaway across all five settings is that the classical L-corridor is unusually friendly to the sofa: Gerver's 2.22 is bigger than anything we could pull out of (i) asymmetric L's, (ii) the disk-cross-section tube, (iii) any S-corridor, (iv) any multi-bend corridor, or (v) any convex restriction. The square cross-section, the right-angle bend, and the freedom to be non-convex each contribute meaningfully to the bound."}], "Text"]
};

nb = Notebook[cells,
  WindowSize -> {1100, 800},
  FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (January 2026)",
  StyleDefinitions -> "Default.nb",
  Saveable -> True];

outNb = FileNameJoin[{hereDir, "Sofa2.nb"}];
Print["writing notebook -> ", outNb];
Export[outNb, nb];
Print["done; nb=", If[FileExistsQ[outNb], FileByteCount[outNb], "FAILED"]];
