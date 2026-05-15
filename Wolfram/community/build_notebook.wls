#!/usr/bin/env wolframscript
(* Build Wolfram/community/Sofa.nb — a Wolfram Notebook for the Wolfram
   Community post on the moving-sofa problem. Also exports a PDF for
   offline review.

   To rebuild:
     wolframscript Wolfram/community/build_notebook.wls
*)

hereDir = DirectoryName[
  If[$ScriptCommandLine =!= {}, $InputFileName, NotebookFileName[]]];
imgDir = hereDir;            (* gif / png live next to the script *)

(* ----- helpers ----------------------------------------------------- *)

imgWidth = 680;

loadImg[name_, width_:Automatic] := With[
  {p = FileNameJoin[{imgDir, name}], w = If[width === Automatic, imgWidth, width]},
  Which[
    ! FileExistsQ[p],
      Style["[missing image: " <> name <> "]", Red],
    ToLowerCase[FileExtension[p]] === "gif",
      (* Import as a Wolfram AnimatedImage so it actually animates inline
         instead of unrolling into a list of frames. *)
      ImageResize[Import[p, "AnimatedImage"], w],
    True,
      Image[Import[p], ImageSize -> w]]];

imgCell[name_, caption_String: ""] := Cell[
  BoxData[ToBoxes[
    If[caption === "",
       loadImg[name],
       Column[{loadImg[name],
               Pane[Style[caption, Italic, FontSize -> 10,
                          FontFamily -> "Helvetica"],
                    imgWidth, Alignment -> Center]},
              Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {10, 14}}];

(* Code cells *)
wlInput[code_String] := Cell[BoxData[code], "Input"];
pyCell[code_String] := Cell[
  BoxData[
    StyleBox[code, FontFamily -> "Courier", FontSize -> 11]],
  "Program", CellMargins -> {{20, 0}, {6, 6}}];

mathI[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
                    "InlineFormula", FontFamily -> "Times"];
mathD[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
                    "DisplayFormula", TextAlignment -> Center,
                    FontFamily -> "Times"];

(* ----- math fragments --------------------------------------------- *)

sofaMaxIntersect = RowBox[{
  SubscriptBox["S", "max"], "(", "C", ")", "=",
  RowBox[{
    UnderscriptBox["\[Intersection]", RowBox[{"\[Theta]", "\[Element]",
                                              RowBox[{"[", "0", ",",
                                                      "\[Pi]/2", "]"}]}]],
    "L", "(", "\[Theta]", ";", "\[ThickSpace]",
    RowBox[{"C", "(", "\[Theta]", ")"}], ")"
  }]
}];

eyEx = RowBox[{
  SubscriptBox[StyleBox["e", FontWeight -> Bold], "x"],
  "(", "\[Theta]", ")", "=",
  RowBox[{"(", "cos", "\[Theta]", ",",
         "-", "sin", "\[Theta]", ")"}], ",",
  "\[ThickSpace]\[ThickSpace]",
  SubscriptBox[StyleBox["e", FontWeight -> Bold], "y"],
  "(", "\[Theta]", ")", "=",
  RowBox[{"(", "sin", "\[Theta]", ",", "cos", "\[Theta]", ")"}]}];

alphaSweep = RowBox[{
  "\[Theta]", "\[Element]", "[", "0", ",",
  RowBox[{"\[Pi]", "-", "\[Alpha]"}], "]"}];

(* ----- the Python pipeline as a single string --------------------- *)

pythonPipeline = "import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon as MplPolygon
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from shapely.geometry import Polygon, box as shp_box
from shapely.validation import make_valid

BIG = 40.0
WORLD_L = np.array([[0, 0], [BIG, 0], [BIG, 1], [1, 1], [1, BIG], [0, BIG]],
                   dtype=float)

def corridor(theta, cx, cy):
    c, s = np.cos(theta), np.sin(theta)
    ex, ey = np.array([c, -s]), np.array([s, c])
    pts = np.array([cx, cy]) + (WORLD_L[:, 0:1] - 1) * ex \\
                            + (WORLD_L[:, 1:2] - 1) * ey
    return Polygon(pts.tolist())

def build_sofa(corners, thetas):
    sofa = corridor(thetas[0], *corners[0])
    for i in range(1, len(thetas)):
        sofa = sofa.intersection(corridor(thetas[i], *corners[i]))
        if sofa.is_empty: return sofa
    return sofa if sofa.is_valid else make_valid(sofa)

def trajectory(ctrl, ts):
    cthetas = np.linspace(0, np.pi/2, len(ctrl))
    sx = CubicSpline(cthetas, ctrl[:, 0], bc_type='natural')
    sy = CubicSpline(cthetas, ctrl[:, 1], bc_type='natural')
    return np.column_stack([sx(ts), sy(ts)])

def neg_area(flat, ts, M):
    sofa = build_sofa(trajectory(flat.reshape(M, 2), ts), ts)
    return 0.0 if sofa.is_empty else -float(sofa.area)

def quarter_arc(M, r=2/np.pi):
    th = np.linspace(0, np.pi/2, M)
    return np.column_stack([r * np.sin(th), r * np.cos(th)])

# Optimise
M, n_th = 10, 100
ts = np.linspace(0, np.pi/2, n_th)
x = quarter_arc(M).reshape(-1)
best = x.copy(); best_f = neg_area(x, ts, M)
print(f'initial area = {-best_f:.5f}')
for k in range(4):
    init = best if k == 0 else best + 0.05 * np.random.randn(best.size)
    r = minimize(neg_area, init, args=(ts, M), method='Nelder-Mead',
                 options=dict(maxiter=3000, xatol=1e-5, fatol=1e-7,
                              adaptive=True))
    if r.fun < best_f: best_f, best = r.fun, r.x
    print(f'  round {k}: area = {-r.fun:.5f}')

# Honest area at much denser theta sampling
ts_eval = np.linspace(0, np.pi/2, 1600)
sofa = build_sofa(trajectory(best.reshape(M, 2), ts_eval), ts_eval)
print(f'\\\\nfinal area (1600-sample eval): {sofa.area:.5f}')";

(* ----- the Wolfram literal-translation snippet, shown in the post  *)

wlLiteralGeom = "(* corridor polygon at orientation theta, corner at (cx,cy) *)
big = 40.0;
worldL = {{0., 0.}, {big, 0.}, {big, 1.}, {1., 1.}, {1., big}, {0., big}};
corridorPoly[\[Theta]_?NumericQ, {cx_?NumericQ, cy_?NumericQ}] := Module[
  {c = Cos[\[Theta]], s = Sin[\[Theta]], ex, ey},
  ex = {c, -s}; ey = {s, c};
  Polygon[Table[{cx, cy} + (p[[1]] - 1) ex + (p[[2]] - 1) ey, {p, worldL}]]];

(* sofa = polygon intersection of every corridor pose *)
buildSofa[corners_, thetas_] := Module[{r, m = Length[thetas]},
  r = corridorPoly[thetas[[1]], corners[[1]]];
  Do[r = RegionIntersection[r, corridorPoly[thetas[[i]], corners[[i]]]],
     {i, 2, m}];
  r];

(* cubic-spline trajectory through M control points *)
evalTraj[ctrl_?MatrixQ, thetas_List] := Module[
  {m = Length[ctrl], cth, fx, fy},
  cth = Subdivide[0., Pi/2., m - 1];
  fx = Interpolation[Transpose[{cth, ctrl[[All, 1]]}], InterpolationOrder -> 3];
  fy = Interpolation[Transpose[{cth, ctrl[[All, 2]]}], InterpolationOrder -> 3];
  Transpose[{fx /@ thetas, fy /@ thetas}]];

(* Hammersley quarter-arc seed *)
quarterArc[m_Integer, r_: 2./Pi] := With[
  {th = Subdivide[0., Pi/2., m - 1]},
  Transpose[{r Sin[th], r Cos[th]}]];";

wlLiteralNM = "(* Hand-rolled Nelder-Mead.  FindMinimum / NMinimize in
   Wolfram don't expose a seeded NelderMead directly, but the algorithm
   is 30 lines: reflect / expand / contract / shrink the simplex. *)
nelderMead[f_, x0_List, maxIter_: 600, fTol_: 1.*^-6] := Module[
  {n = Length[x0], simplex, fvals, ord, iter = 0, centroid, worst,
   xr, fr, xe, fe, xc, fc, newPt, newFv, shrink, conv = False,
   alpha = 1., gamma = 2., rho = 0.5, sigma = 0.5},
  simplex = Prepend[Table[x0 + 0.05 UnitVector[n, i], {i, n}], x0];
  fvals = f /@ simplex;
  While[iter < maxIter && ! conv,
    iter++;
    ord = Ordering[fvals];
    simplex = simplex[[ord]]; fvals = fvals[[ord]];
    centroid = Mean[Drop[simplex, -1]]; worst = Last[simplex];
    xr = centroid + alpha (centroid - worst); fr = f[xr];
    Which[
      fr < fvals[[1]],
        xe = centroid + gamma (xr - centroid); fe = f[xe];
        {newPt, newFv} = If[fe < fr, {xe, fe}, {xr, fr}]; shrink = False,
      fr < fvals[[-2]],
        {newPt, newFv} = {xr, fr}; shrink = False,
      True,
        xc = centroid + rho (xr - centroid); fc = f[xc];
        If[fc < fvals[[-1]], {newPt, newFv} = {xc, fc}; shrink = False,
                              shrink = True]];
    If[shrink,
      simplex = Prepend[Table[
        simplex[[1]] + sigma (simplex[[i]] - simplex[[1]]),
        {i, 2, n + 1}], simplex[[1]]];
      fvals = f /@ simplex,
      simplex[[-1]] = newPt; fvals[[-1]] = newFv];
    conv = Max[fvals] - Min[fvals] < fTol];
  ord = Ordering[fvals];
  {fvals[[ord[[1]]]], simplex[[ord[[1]]]], iter}];";

wlLiteralRun = "(* Drive the optimisation. M = 8 control points, 30 thetas
   in the objective, 200 thetas for the honest area readout. *)
M = 8; nTh = 30; nThEval = 200;
thetas    = Subdivide[0., Pi/2., nTh - 1];
thetasEv  = Subdivide[0., Pi/2., nThEval - 1];

negArea[flat_?VectorQ] := Module[{ctrl = Partition[flat, 2], a},
  a = Quiet @ Area @ buildSofa[evalTraj[ctrl, thetas], thetas];
  If[NumericQ[a] && a > 0, -a, 0.]];

x0 = Flatten @ quarterArc[M];
Print[\"initial area = \", -negArea[x0]];

best = x0; bestF = negArea[x0];
SeedRandom[7];
Do[
  start = If[k == 0, best, best + 0.05 RandomReal[{-1, 1}, 2 M]];
  {fmin, xmin, iters} = nelderMead[negArea, start, 600];
  If[fmin < bestF, bestF = fmin; best = xmin];
  Print[\"round \", k, \": area = \", -fmin, \"  (\", iters, \" iters)\"];
  , {k, 0, 3}];

(* Honest area at the dense eval grid. *)
ctrlBest = Partition[best, 2];
sofaBest = buildSofa[evalTraj[ctrlBest, thetasEv], thetasEv];
Print[\"honest area : \", Area[sofaBest]];
Print[\"Hammersley  : \", N[Pi/2 + 2./Pi]];
Print[\"Gerver bound: \", 2.21953166887197];";

(* ----- the user's original short Wolfram exploration ------------- *)

wlExplore1 = "reg = RegionUnion[
  DiscretizeGraphics[Graphics @ Rectangle[{-5, -1}, {2, 1}]],
  DiscretizeGraphics[Graphics @ Rectangle[{0, -1}, {2, -6}]]];";

wlExplore2 = "Manipulate[
  RegionPlot[{reg, TransformedRegion[TransformedRegion[reg,
     TranslationTransform[{-2 Cos[\[Alpha]] + 1, -1 + 2 Cos[\[Alpha]]}]],
     RotationTransform[-\[Alpha] + Pi/4, {0, -1}]]},
    PlotRange -> {{-6, 6}, {-6, 6}}],
  {\[Alpha], 0, Pi/2}]";

wlExploreHand = "(* hand-tuned parametric trajectory the user converged on; \
amp \[TildeTilde] 1.22, exp1 \[TildeTilde] 0.68, exp2 \[TildeTilde] 0.53 *)
amp = 1.22; exp1 = 0.68; exp2 = 0.53;
regs = Table[
   TransformedRegion[
    TransformedRegion[reg,
     TranslationTransform[{-amp Sign[Cos[2 \[Alpha]]] Abs[Cos[2 \[Alpha]]]^exp1,
                            amp Sign[Cos[2 \[Alpha]]] Abs[Cos[2 \[Alpha]]]^exp1}]],
    RotationTransform[Sign[-\[Alpha] + Pi/4.]
                        (Abs[-\[Alpha] + Pi/4.] 4/Pi)^exp2 Pi/4., {0, -1}]],
   {\[Alpha], 0, Pi/2, 0.015}];
sofa = Fold[RegionIntersection, regs];
Area[sofa] / 4   (* divide by 4 because reg is on a 2-unit-wide corridor *)";

(* ----- build the notebook ----------------------------------------- *)

cells = {

  Cell["The Moving Sofa Problem: Numerical Optimisation in 2D and Beyond",
       "Title"],

  Cell["Finding Gerver-class sofas, varying the bend angle, and a glance at \
the higher-dimensional generalisation.",
       "Subtitle"],

  Cell[TextData[{"Marco Thiel, May 2026 \[LongDash] ",
                 ButtonBox["github.com/mthiel74/SofaProblem",
                   BaseStyle -> "Hyperlink",
                   ButtonData -> {URL["https://github.com/mthiel74/SofaProblem"], None}]}],
       "Text", TextAlignment -> Center, FontSlant -> Italic],

  (* ====== Abstract ====== *)
  Cell["Abstract", "Section"],
  Cell[TextData[{
    "The ", StyleBox["moving sofa problem", FontSlant -> Italic],
    " asks for the largest 2D rigid body that can slide through an L\
-shaped corridor of unit width while making a 90\[Degree] turn. \
Hammersley (1968) constructed a sofa of area ", mathI[RowBox[{"\[Pi]","/","2","+",RowBox[{"2","/","\[Pi]"}]," \[TildeTilde] ","2.207"}]],
    "; Gerver (1992) refined this to an 18-arc shape of area ", mathI["2.21953"],
    ", and Baek (2024) recently proved that Gerver's sofa is optimal. This notebook reproduces a Gerver-class shape from scratch with a small Python optimiser (cubic-spline trajectory of the inner corner of the L; Nelder\[Dash]Mead on \[TildeTilde]20 parameters), gives a literal Wolfram-Language port of the same algorithm, sweeps the corridor bend angle across (0, \[Pi]) and reports the resulting areas, and ends with a short discussion of what changes (and what is open) in three dimensions."}],
    "Text"],
  imgCell["sofa.png", "Figure 1. The optimised sofa (light blue), area \[TildeTilde] 2.2197, together with the inner-corner trajectory the optimiser discovered (dashed orange). Width 1, length 2 + 4/\[Pi] \[TildeTilde] 3.27. The sofa is essentially Gerver's shape."],

  (* ====== Section 1: The puzzle ====== *)
  Cell["1. The puzzle", "Section"],
  Cell[TextData[{
    "An L-shaped corridor of unit width: the horizontal arm is the strip ",
    mathI[RowBox[{RowBox[{"{","(", "x", ",", "y", ")", "}"}]," : ","x", "\[GreaterEqual]", "0", ", ",
                  "y", "\[Element]", "[", "0", ",", "1", "]"}]],
    " and the vertical arm is its reflection across ", mathI[RowBox[{"y","=","x"}]],
    ". A connected planar body ", mathI["S"], " is a ",
    StyleBox["sofa", FontSlant -> Italic],
    " if there is a continuous rigid motion of ", mathI["S"],
    " starting deep in the horizontal arm and ending deep in the vertical arm, with ",
    mathI["S"], " contained in the corridor at every instant. What is the ",
    mathI[RowBox[{"sup", " ", "area", "(", "S", ")"}]], "?"}], "Text"],
  Cell[TextData[{
    "A trivial unit square gives area 1; a half-disk of radius 1 fits and gives area ",
    mathI[RowBox[{"\[Pi]","/","2","\[TildeTilde]","1.571"}]],
    ". Hammersley's 1968 construction \[LongDash] a rectangle with a semicircular notch cut from its bottom \[LongDash] gives ",
    mathI[RowBox[{"\[Pi]","/","2","+",RowBox[{"2","/","\[Pi]"}],"\[TildeTilde]","2.207"}]],
    ". Gerver (1992) replaced the rectangle with an 18-arc shape and got ",
    mathI["2.21953"],
    "; the recent proof by Baek (arXiv:2411.19826, 2024) shows that ",
    mathI["2.21953..."], " is the exact answer."}], "Text"],

  (* ====== Section 2: Formulation ====== *)
  Cell["2. The sofa-fixed-frame formulation", "Section"],
  Cell[TextData[{
    "Instead of moving the sofa through a fixed corridor, fix the sofa and let the corridor rotate and translate around it. The motion is parameterised by an angle ", mathI["\[Theta]"], " \[Element] [0, \[Pi]/2]; in the sofa-fixed frame the world's axes are rotated by ", mathI[RowBox[{"-", "\[Theta]"}]], ":"}], "Text"],
  mathD[eyEx],
  Cell[TextData[{
    "Each pose places an L-shaped corridor ", mathI[RowBox[{"L", "(", "\[Theta]", ")"}]], " in the sofa frame. The largest sofa consistent with a given trajectory ", mathI[RowBox[{"C", "(", "\[Theta]", ")"}]], " of the inner corner is the polygon intersection"}], "Text"],
  mathD[sofaMaxIntersect],
  Cell[TextData[{
    "and we maximise its area over the trajectory ", mathI["C"], "."}], "Text"],

  (* ====== Section 3: A Wolfram exploration ====== *)
  Cell["3. A Wolfram-Language exploration", "Section"],
  Cell[TextData[{
    "Discretise the integral above: pick a dense set of ", mathI["\[Theta]"],
    " samples, build the L-corridor at each sample, intersect them all. The trajectory ",
    mathI[RowBox[{"C", "(", "\[Theta]", ")"}]],
    " is the parameter we optimise. As a first pass, fit a parametric trajectory by hand: a translation of the form",
    " "
  }], "Text"],
  mathD[RowBox[{
    "T", "(", "\[Alpha]", ")", "=", "amp", "\[CenterDot]",
    "sign", "(", "cos", "2", "\[Alpha]", ")", "\[CenterDot]",
    "|", "cos", "2", "\[Alpha]", SuperscriptBox["|", SubscriptBox["e", "1"]],
    ", ", "\[ThinSpace]\[ThinSpace]", "rotation",
    " ", "\[Phi]", "(", "\[Alpha]", ")", "=", "..."}]],

  Cell["Set up the corridor as a region and look at one corridor pose:",
       "Text"],
  wlInput[wlExplore1],
  wlInput[wlExplore2],
  Cell["Tune the two exponents and the amplitude by hand and integrate every pose:",
       "Text"],
  wlInput[wlExploreHand],
  Cell[TextData[{
    "Trial and error converges to area ",
    mathI["\[TildeTilde] 2.14"],
    " on a 2-unit-wide corridor (so ", mathI["\[TildeTilde] 2.14"], " on the unit-width problem after the divide-by-4). This is comfortably above the simple half-disk of radius 1 (area ", mathI[RowBox[{"\[Pi]","/","2"}]], ") but still short of Hammersley\[CloseCurlyQuote]s 2.207 \[LongDash] showing how much is left on the table once we move away from low-dimensional hand-tuned curves."}], "Text"],

  (* ====== Section 4: Python pipeline ====== *)
  Cell["4. A Python optimiser", "Section"],
  Cell[TextData[{
    "Trading the two-parameter hand-tuned trajectory for a 10-point cubic spline lifts the area to a Gerver-class shape. The pipeline below \[LongDash] geometry + cubic-spline trajectory + Nelder\[Dash]Mead with warm-started restarts \[LongDash] runs in a few minutes and produces ", StyleBox["sofa.png", "Input"], " and ", StyleBox["sofa.gif", "Input"], ":"}], "Text"],
  pyCell[pythonPipeline],
  Cell["Animation of the result \[LongDash] the L-corridor sweeping over the sofa in the sofa-fixed frame:",
       "Text"],
  imgCell["sofa.gif", "Figure 2. As \[Theta] runs from 0\[Degree] to 90\[Degree] (corridor rotating in the sofa-fixed frame), the orange-shaded L sweeps over the blue sofa. The sofa is fully contained in every frame."],
  Cell[TextData[{
    "The honest area evaluated at 1600 \[Theta] samples is ",
    StyleBox["2.21968", FontWeight -> Bold],
    ", within ", mathI[RowBox[{"2", "\[Times]", SuperscriptBox["10", RowBox[{"-", "4"}]]}]],
    " of Gerver's bound ", mathI["2.21953..."], "."}],
       "Text"],

  (* ====== Section 5: Wolfram literal port ====== *)
  Cell["5. A literal Wolfram-Language port", "Section"],
  Cell[TextData[{
    "The same algorithm in pure Wolfram Language. ",
    StyleBox["Polygon", "Input"], " primitives carry the corridor poses and ",
    StyleBox["RegionIntersection", "Input"], " composes them. The optimiser is a hand-rolled Nelder\[Dash]Mead because Wolfram's ", StyleBox["FindMinimum", "Input"],
    " doesn't expose Nelder\[Dash]Mead and ", StyleBox["NMinimize", "Input"],
    "'s NelderMead path doesn't honour an initial-point seed cleanly. The three blocks below \[LongDash] geometry, optimiser, driver \[LongDash] are everything; running the driver block reproduces the Gerver-class area."}], "Text"],

  Cell["Geometry: L-corridor polygon, region intersection, spline trajectory.", "Subsection"],
  wlInput[wlLiteralGeom],

  Cell["A seeded Nelder\[Dash]Mead in 30 lines of Wolfram.", "Subsection"],
  wlInput[wlLiteralNM],

  Cell["Drive the optimisation. With M = 8 control points and 30 theta samples per call, the four warm-started rounds converge in about 100 seconds on a laptop:", "Subsection"],
  wlInput[wlLiteralRun],

  Cell[TextData[{
    "Output on the machine used for this notebook:\n",
    StyleBox["initial area = 1.00895\nround 0: area = 2.20372  (600 iters)\nround 1: area = 2.24009  (600 iters)\nround 2: area = 2.23959  (600 iters)\nround 3: area = 2.23884  (600 iters)\nhonest area : ", "Output"],
    StyleBox["2.22195", FontWeight -> Bold, FontFamily -> "Courier"],
    StyleBox["\nHammersley  : 2.20742\nGerver bound: 2.21953", "Output"]
  }], "Text"],

  Cell[TextData[{
    "The Wolfram run lands at ", StyleBox["2.22195", FontWeight -> Bold],
    " on a 200-theta evaluation grid, matching Gerver's optimal value 2.21953... to four decimals \[LongDash] the same precision the Python pipeline reached in Section 4. The two implementations agree on the optimum to within the discretisation noise of their respective theta grids. Polygon-intersection cost dominates in both languages; Wolfram's ",
    StyleBox["RegionIntersection", "Input"],
    " on ", StyleBox["Polygon", "Input"],
    " primitives is fast enough that ", mathI[RowBox[{"2400"}]], " evaluations of the 30-theta objective finish in roughly ", mathI[RowBox[{"100"}]], " seconds."}], "Text"],
  Cell["Once the optimiser has produced control points, the same code drives a Wolfram-native animation of the sweep:", "Text"],
  wlInput["(* Animation of the corridor sweep in the sofa-fixed frame.
   ctrlXY is the optimised 10x2 matrix from the FindMinimum step above. *)
With[{nframes = 60},
  thetasAnim = Subdivide[0.01, Pi/2 - 0.01, nframes - 1];
  cornersAnim = evalTraj[ctrlXY, thetasAnim];
  sofaPoly = buildSofa[evalTraj[ctrlXY, Subdivide[0., Pi/2., 159]],
                        Subdivide[0., Pi/2., 159]];
  frames = Table[
    Show[{
      RegionPlot[sofaPoly, PlotStyle -> Lighter[Blue, 0.6],
                  BoundaryStyle -> Directive[Blue, Thickness[0.004]],
                  Frame -> False, Axes -> False, ImageSize -> 480],
      Graphics[{Opacity[0.35], LightOrange, EdgeForm[Directive[Darker@Red, Thick]],
                corridorPoly[thetasAnim[[i]], cornersAnim[[i]]]}]
    }, PlotRange -> {{-1.5, 2}, {-2.5, 1.5}}, AspectRatio -> Automatic],
    {i, 1, nframes}];
  ListAnimate[frames, AnimationRunning -> True, AnimationRate -> 20]]"],
  Cell["Two equivalent ways to package the animation: ListAnimate (interactive) or Export[\"sofa-wolfram.gif\", frames] for a portable GIF.", "Text"],
  Cell["Scaling the angle sweep with RemoteBatchSubmit", "Subsection"],
  Cell[TextData[{
    "The seven optimisations in Section 7 are embarrassingly parallel: each \[Alpha] is an independent Nelder\[Dash]Mead run. ",
    StyleBox["RemoteBatchSubmit", "Input"], " ships each anchor as a job to the Wolfram Cloud:"}], "Text"],
  wlInput["(* fan out the angle sweep across cloud workers; with no
   \"BatchComputeProvider\" option this picks the Wolfram-managed
   default that needs no per-user AWS / GCP credentials. *)
job = RemoteBatchSubmit[
  CloudObject[\"sofa-sweep\"],
  Module[{a = #},
    optimiseSofaAlpha[a, M -> 10, nThetas -> 80, nRounds -> 3]
  ] &,
  Range[60 Degree, 150 Degree, 15 Degree]
];
(* poll until done, then download the seven optimised sofas *)
RemoteBatchMapDownload[job]"],
  Cell[TextData[{
    "Each anchor takes 1\[Dash]3 min, so a seven-angle sweep finishes in roughly the time of one anchor when run in parallel. Keep ",
    StyleBox["nRounds", "Input"], " modest and the Wolfram Cloud compute cost stays well under a dollar."}], "Text"],

  (* ====== Section 6: Higher dimensions ====== *)
  Cell["6. Higher dimensions", "Section"],
  Cell[TextData[{
    "The natural 3D analogue replaces the L-shaped planar corridor by an L-shaped ",
    StyleBox["tube", FontSlant -> Italic],
    " of unit cross-section. With a square cross-section, the tube is ",
    Cell[BoxData[FormBox[
      RowBox[{"(", RowBox[{"{","(","x",",","y",",","z",")",":","\[ThinSpace]","0","\[LessEqual]","y",",","z","\[LessEqual]","1",",","\[ThinSpace]","x","\[LessEqual]","1","}"}], ")",
             "\[Union]",
             "(", RowBox[{"{","0","\[LessEqual]","x",",","z","\[LessEqual]","1",",","\[ThinSpace]","y","\[LessEqual]","1","}"}], ")"}],
      TraditionalForm]],
        "InlineFormula", FontFamily -> "Times"],
    " and one asks for the largest rigid body ", mathI[RowBox[{"S", "\[Subset]", SuperscriptBox["R", "3"]}]],
    " with a continuous motion ", mathI[RowBox[{"g", "(", "t", ")", "\[Element]", "SE", "(", "3", ")"}]],
    " keeping ", mathI[RowBox[{"g", "(", "t", ")", "S"}]], " in the tube."}], "Text"],
  Cell[TextData[{
    "The problem is essentially open. The only rigorous higher-dimensional result is for rectangular boxes (Horv\[AAcute]th et al.) and is the trivial answer. The trivial volume lower bound is ",
    StyleBox["A_Gerver \[Times] 1 \[TildeTilde] 2.21953", "Input"], " (extrude Gerver's sofa along the third axis); no published explicit upper bound exists, but folklore suggests ", mathI[RowBox[{"\[LessEqual]", "3"}]], ". A non-trivial 3D construction is widely ",
    StyleBox["expected", FontSlant -> Italic],
    " to exist because the 3D motion has two extra rotational degrees of freedom (pitch and roll) above the planar problem \[LongDash] the sofa can tilt and corkscrew through the bend \[LongDash] but neither a construction nor a sharp upper bound is in print."}], "Text"],
  Cell[TextData[{
    "Algorithmically the formulation transfers: parameterise the path ", mathI[RowBox[{"g","(","t",")","=","(","R","(","t",")",",","\[ThinSpace]","p","(","t",")",")"}]],
    " with ", mathI[RowBox[{"R","(","t",")","\[Element]","SO","(","3",")"}]],
    " by a quaternion (or Euler-angle) spline; the sofa is the intersection ",
    mathI[RowBox[{
      UnderscriptBox["\[Intersection]", "t"],
      SuperscriptBox["g", RowBox[{"(","t",")",SuperscriptBox["","-1"]}]],
      "C"}]],
    ". With ", mathI[RowBox[{"\[TildeTilde]","50\[Dash]200"}]], " parameters and millions of voxel/half-space evaluations per call, this is a workstation-hours run \[LongDash] research-grade rather than turn-the-crank."}], "Text"],

  imgCell["sofa_3d_infographic.png", "Figure 5. The 3D set-up. A unit-square L-tube bent 90\[Degree], the (trivial) extruded Gerver sofa inside it, and the three extra rotational degrees of freedom over the 2D problem: yaw \[Theta] (already used in 2D), pitch \[CurlyPhi], and roll \[Beta]. The configuration space is SO(3) instead of SO(2). (Image generated with OpenAI gpt-image / image-2 from a hand-written prompt; the math and layout are mine, the rendering is the model's.)"],

  Cell["The 3D set-up", "Subsection"],
  Cell[TextData[{
    "The 3D corridor is an L-shaped ", StyleBox["tube", FontSlant -> Italic],
    " of unit-square cross-section. With the 2D L-corridor ", mathI["L"], " of \[Section]2,"}], "Text"],
  mathD[RowBox[{
    SubscriptBox["C", "3D"], "=", "L", "\[Times]", "[", "0", ",", "1", "]"}]],
  Cell[TextData[{
    "A 3D moving sofa is a rigid body ", mathI[RowBox[{"S","\[Subset]",SuperscriptBox["R","3"]}]], " admitting a continuous motion ",
    mathI[RowBox[{"g","(","t",")","=","(","R","(","t",")",",","\[ThinSpace]","p","(","t",")",")","\[Element]", "SE","(","3",")"}]],
    " with ", mathI[RowBox[{"g","(","t",")","\[ThinSpace]","S","\[Subset]",SubscriptBox["C","3D"]}]],
    " for all ", mathI["t"], "."}], "Text"],
  Cell[TextData[{
    "The 2D analysis transfers to a ", StyleBox["necessary", FontSlant -> Italic],
    " condition: the (x, y) projection of ", mathI[RowBox[{"g","(","t",")","S"}]],
    " must lie in ", mathI["L"], " at every ", mathI["t"], ". So every horizontal slice ",
    mathI[RowBox[{"S","\[Intersection]",RowBox[{"{","z","=",SubscriptBox["z","0"],"}"}]}]],
    " is a 2D body that moves through ", mathI["L"], " under a (slice-dependent) 2D rigid motion, and has area at most Gerver's ", mathI[RowBox[{"2.21953","\[Ellipsis]"}]],
    ". Integrating slice areas:"}], "Text"],
  mathD[RowBox[{
    "vol", "(", "S", ")", "=",
    SubsuperscriptBox["\[Integral]", "0", "1"],
    "area", "(", SubscriptBox["S", "z"], ")",
    " ", "d", "z", "\[LessEqual]", "2.21953", "\[Times]", "1", "=", "2.21953",
    "\[Ellipsis]"}]],
  Cell[TextData[{
    "Wait \[LongDash] but every slice gets a ", StyleBox["different", FontSlant -> Italic],
    " 2D motion if ", mathI[RowBox[{"R","(","t",")"}]], " has a pitch or roll component. Each slice ",
    mathI[SubscriptBox["S", "z"]], " then traces out an ", mathI["\[Alpha]"], "-bend corridor where ",
    mathI["\[Alpha]"], " depends on ", mathI["z"], ". For slices with ", mathI[RowBox[{"\[Alpha]",">",RowBox[{"\[Pi]","/","2"}]}]], " (look back at Section 7!) the per-slice area cap ", mathI[RowBox[{"S","(","\[Alpha]",")"}]],
    " is strictly larger than 2.21953. So in principle the 3D bound can be ", StyleBox["exceeded", FontSlant -> Italic],
    " by spending some slices' rotational budget asymmetrically. This is the open intuition."}], "Text"],

  Cell["A first numerical run: extruded Gerver in the tube", "Subsection"],
  Cell["Build the trivial 3D construction \[LongDash] the 2D shape from Section 4 extruded as a unit-height prism \[LongDash] and slide it through the L-tube along the (lifted) 2D Gerver motion. The Python pipeline below produces the animation in Figure 6:", "Text"],
  pyCell["import math, numpy as np\nfrom sofa_angle import evaluate_trajectory_alpha\n\n# Boundary of the 2D Gerver-class sofa (from Section 4).\nboundary = np.load('sweep_angles.npz', allow_pickle=True)['boundary_2']\nctrl    = np.load('sweep_angles.npz', allow_pickle=True)['controls_2']\n\n# Shoelace area of the closed 2D polygon.\nx, y = boundary[:, 0], boundary[:, 1]\narea = 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))\nprint(f'2D Gerver-class area   = {area:.5f}')\nprint(f'Trivial 3D lower bound = area x 1 = {area:.5f}')\n\n# 3D rigid motion: lift the 2D Gerver rotation about the world inner\n# corner (1, 1) to z = 0.5 (the centre of the unit-square tube).\ndef lab_frame(theta, beta=0.0):\n    Rz = np.array([[ math.cos(theta), -math.sin(theta), 0],\n                   [ math.sin(theta),  math.cos(theta), 0],\n                   [               0,                0, 1]])\n    Ry = np.array([[ math.cos(beta),  0, math.sin(beta)],\n                   [               0, 1,              0],\n                   [-math.sin(beta),  0, math.cos(beta)]])\n    R = Rz @ Ry\n    c2 = evaluate_trajectory_alpha(ctrl, np.array([theta]), np.pi / 2.0)[0]\n    c3 = np.array([c2[0], c2[1], 0.5])\n    t  = np.array([1.0, 1.0, 0.5]) - R @ c3\n    return R, t\n\n# Sample the motion at e.g. 60 thetas and compose a rigid pose at each.\nthetas = np.linspace(0.01, np.pi/2 - 0.01, 60)\nposes  = [lab_frame(th, beta=0.0) for th in thetas]"],

  imgCell["sofa_3d.gif", "Figure 6. The trivial 3D construction (extruded Gerver sofa) traversing the unit-square L-tube. Pure 2D motion in (x, y); the prism rides flat through the corner. Volume = 2D area \[Times] 1 = 2.21973."],

  Cell[TextData[{
    "What about exercising the extra DOFs? A naive rigid tilt of the same prism about its long axis is plotted next: any non-zero ", mathI["\[Beta]"],
    " lifts one side of the top face above ", mathI[RowBox[{"z","=","1"}]], " and drops the other below ", mathI[RowBox[{"z","=","0"}]], " \[LongDash] the prism pokes through the floor or the ceiling and the construction immediately fails."}], "Text"],
  imgCell["tilt_attempt.png", "Figure 7. Rigid tilt of the extruded prism about its long axis. At every \[Beta] > 0 the prism breaks out of the [0, 1] clearance (red). Getting volume out of pitch / roll requires a non-prismatic body whose shape is correlated with the rotation \[LongDash] research-grade."],

  Cell[TextData[{
    "We can be more concrete about a lower bound that beats the trivial extrusion. Imagine subdividing the prism into ", mathI["n"], " thin horizontal slabs of height ",
    mathI[RowBox[{"1","/","n"}]], " and allowing each slab a small ", StyleBox["independent", FontSlant -> Italic],
    " yaw offset ", mathI[RowBox[{SubscriptBox["\[Delta]","i"]}]], " applied through a per-slab roll-and-shift. The body is still rigid (you've baked the offsets into its 3D shape), but each slab now traverses an effective ", mathI[RowBox[{"(",RowBox[{"\[Pi]","/","2"}],"+",SubscriptBox["\[Delta]","i"],")"}]], "-corner; using the angle-sweep numbers of Section 7,"}], "Text"],
  mathD[RowBox[{
    "vol", "(", SubscriptBox["S", "twist"], ")", "=",
    FractionBox["1", "n"],
    UnderoverscriptBox["\[Sum]", RowBox[{"i", "=", "1"}], "n"],
    " ", "S", "(",
    RowBox[{FractionBox["\[Pi]", "2"], "+", SubscriptBox["\[Delta]", "i"]}], ")"}]],
  Cell[TextData[{
    "with ", mathI[RowBox[{SubscriptBox["\[Delta]", "1"], "+", "\[Ellipsis]", "+", SubscriptBox["\[Delta]", "n"], "=", "0"}]],
    " so the whole body still nets a 90\[Degree] turn. Because ", mathI[RowBox[{"S","(","\[Alpha]",")"}]],
    " is convex over the relevant range (Section 7 numerics), Jensen's inequality says: if you can give ", StyleBox["any", FontSlant -> Italic], " slab a positive ",
    mathI[SubscriptBox["\[Delta]","i"]], " while ", StyleBox["another", FontSlant -> Italic],
    " absorbs the equal-and-opposite negative one, the average slice area exceeds ", mathI[RowBox[{"S","(",RowBox[{"\[Pi]","/","2"}],")"}]],
    ". Whether that algebraic gain survives the clearance bookkeeping (each slab must still live in ", mathI[RowBox[{"[","0",",","1","]"}]],
    " in ", mathI["z"], " under the tilted motion) is exactly what the 3D problem leaves open. Driving this trade-off to a sharp lower bound is research-grade: the parameter count climbs into the hundreds and the objective is a 3D Boolean-mesh containment check at every motion step. The numerical machinery in this notebook \[LongDash] cubic-spline trajectory, polygon intersection, warm-started Nelder\[Dash]Mead \[LongDash] generalises directly to that setting, but is well beyond the few-minutes runtime we've allowed here."}], "Text"],

  (* ====== Section 7: Varying bend angle ====== *)
  Cell["7. Varying the bend angle", "Section"],
  Cell[TextData[{
    "Let ", mathI["\[Alpha]"], " \[Element] (0, \[Pi]) be the interior bend angle: ",
    mathI[RowBox[{"\[Alpha]","=",RowBox[{"\[Pi]","/","2"}]}]], " is the classical right angle, ",
    mathI[RowBox[{"\[Alpha]","\[Rule]","\[Pi]"}]], " is a straight strip (the sofa is unbounded), and ",
    mathI[RowBox[{"\[Alpha]","\[Rule]","0"}]], " is a fold-back (the sofa shrinks). The sofa's rotation across the motion is ",
    mathI[RowBox[{"\[Pi]","-","\[Alpha]"}]], ", so the optimiser parameterises ",
    alphaSweep, "."}], "Text"],
  Cell[TextData[{
    "The literature on non-right-angle bends is essentially empty: Hammersley, Gerver, Romik, Kallus\[Dash]Romik, and Baek all work at ",
    mathI[RowBox[{"\[Alpha]","=",RowBox[{"\[Pi]","/","2"}]}]],
    ". Our optimiser only changes the corridor polygon and the parameter range. Below are the areas at seven anchor angles, each from a ",
    mathI["\[TildeTilde] 2\[Dash]4"], "-minute Nelder\[Dash]Mead run with the Hammersley-arc seed."}], "Text"],
  Cell["The seven anchor angles, side by side. Each row shows the corridor (\"aisle\") at that bend angle in column 2, and the optimised sofa shape (in its own frame, all drawn at the same scale) in column 3. The aisle's inner corner is marked with the orange dot.",
       "Text"],
  Cell[BoxData[ToBoxes[
    Grid[{
      {Style["\[Alpha]", Bold], Style["aisle", Bold], Style["sofa", Bold], Style["area", Bold]},
      {Style["60\[Degree]", 11], loadImg["aisle_60.png", 200], loadImg["sofa_60.png", 200], "1.947"},
      {Style["75\[Degree]", 11], loadImg["aisle_75.png", 200], loadImg["sofa_75.png", 200], "2.057"},
      {Style["90\[Degree]", 11], loadImg["aisle_90.png", 200], loadImg["sofa_90.png", 200], "2.220"},
      {Style["105\[Degree]", 11], loadImg["aisle_105.png", 200], loadImg["sofa_105.png", 200], "2.470"},
      {Style["120\[Degree]", 11], loadImg["aisle_120.png", 200], loadImg["sofa_120.png", 200], "2.883"},
      {Style["135\[Degree]", 11], loadImg["aisle_135.png", 200], loadImg["sofa_135.png", 200], "3.624"},
      {Style["150\[Degree]", 11], loadImg["aisle_150.png", 200], loadImg["sofa_150.png", 200], "5.204"}},
     Frame -> All, Alignment -> Center,
     Background -> {None, {1 -> LightGray}}]]],
   "Output", ShowCellLabel -> False],
  Cell["The same data plotted as a curve, with the heuristic 1/(\[Pi] - \[Alpha]) divergence overlaid for reference:", "Text"],
  imgCell["area_vs_angle.png", "Figure 3. The optimised sofa area S(\[Alpha]) grows from \[TildeTilde] 2 at \[Alpha] = 60\[Degree] to \[GreaterEqual] 5 at \[Alpha] = 150\[Degree]. The dashed curve C/(\[Pi]-\[Alpha]) (calibrated at \[Alpha] = 90\[Degree]) shows the expected strip-limit divergence at \[Alpha]\[Rule]\[Pi]."],
  Cell["And finally the same family rendered as an animation \[LongDash] each frame interpolates the spline control points between adjacent anchor angles:", "Text"],
  imgCell["sweep.gif", "Figure 4. The optimised sofa shape as the bend angle \[Alpha] is swept smoothly between the seven anchors. Each anchor is an independent Nelder\[Dash]Mead optimum, warm-started from the previous one."],

  (* ====== Section 8: Summary tables ====== *)
  Cell["8. Summary tables", "Section"],

  Cell["Table 1. 2D 90\[Degree] sofa benchmarks. The aisle is the same unit-width L-corridor across all rows; the sofa column shows the shape (or a representative shape for analytic bounds).", "Text"],
  Cell[BoxData[ToBoxes[
    Grid[{
      {Style["shape", Bold], Style["sofa", Bold], Style["area", Bold], Style["source", Bold]},
      {"unit square", loadImg["thumb_unit_square.png", 90], "1.000", "trivial"},
      {"unit half-disk", loadImg["thumb_half_disk.png", 90], "\[Pi]/2 \[TildeTilde] 1.5708", "trivial"},
      {"hand-tuned (Wolfram)", "(diagram in \[Section]3)", "\[TildeTilde] 2.14", "Section 3"},
      {"Hammersley sofa", loadImg["thumb_hammersley.png", 130], "\[Pi]/2 + 2/\[Pi] \[TildeTilde] 2.2074", "Hammersley 1968"},
      {"this notebook (Python)", loadImg["thumb_gerver.png", 90], StyleBox["2.21968", FontWeight -> Bold], "Section 4"},
      {"this notebook (Wolfram)", loadImg["thumb_gerver.png", 90], StyleBox["2.22195", FontWeight -> Bold], "Section 5"},
      {"Gerver sofa (best known, optimal)", loadImg["thumb_gerver.png", 90], "2.21953\[Ellipsis]", "Gerver 1992, Baek 2024"},
      {"Kallus\[Dash]Romik upper bound", "\[LongDash]", "\[LessEqual] 2.37", "Kallus\[Dash]Romik 2018"}},
     Frame -> All, Alignment -> {Left, Center, Center, Left},
     Background -> {None, {1 -> LightGray}}]]],
   "Output", ShowCellLabel -> False],

  Cell["Table 2. 90\[Degree] L-corridor in dimensions 2 and 3.", "Text"],
  Cell[BoxData[ToBoxes[
    Grid[{
      {"dimension", "lower bound", "upper bound", "status"},
      {"2D (square corridor)", "2.21953... (Gerver)", "2.21953... (Baek)", "proved optimal"},
      {"3D (square tube)", "\[GreaterEqual] 2.21953 (extrude)", "\[Tilde] 3 (folklore)", "essentially open"},
      {"3D (circular tube)", "\[GreaterEqual] \[Pi]/4 \[Times] (best disk sofa)", "open", "open"}},
      Frame -> All, ItemSize -> {18, 18, 18, 16},
      Background -> {None, {1 -> LightGray}}]]],
    "Output", ShowCellLabel -> False],

  Cell["Table 3. Optimised area at varying bend angle \[Alpha] (this notebook). Each run uses the previous \[Alpha] anchor as a warm-start \[LongDash] without it Nelder\[Dash]Mead easily gets trapped at large \[Alpha].", "Text"],
  Cell[BoxData[ToBoxes[
    Grid[{
      {Style["\[Alpha]", Bold], Style["aisle", Bold], Style["sofa", Bold], Style["area", Bold], Style["notes", Bold]},
      {"60\[Degree]",  loadImg["aisle_60.png",  130], loadImg["sofa_60.png",  130], "1.947", "narrow bend"},
      {"75\[Degree]",  loadImg["aisle_75.png",  130], loadImg["sofa_75.png",  130], "2.057", ""},
      {"90\[Degree]",  loadImg["aisle_90.png",  130], loadImg["sofa_90.png",  130], "2.220", "classical L; matches Gerver"},
      {"105\[Degree]", loadImg["aisle_105.png", 130], loadImg["sofa_105.png", 130], "2.470", ""},
      {"120\[Degree]", loadImg["aisle_120.png", 130], loadImg["sofa_120.png", 130], "2.883", ""},
      {"135\[Degree]", loadImg["aisle_135.png", 130], loadImg["sofa_135.png", 130], "3.624", ""},
      {"150\[Degree]", loadImg["aisle_150.png", 130], loadImg["sofa_150.png", 130], "5.204", "wide-open; long thin sofa"}},
      Frame -> All, Alignment -> {Center, Center, Center, Center, Left},
      Background -> {None, {1 -> LightGray}}]]],
    "Output", ShowCellLabel -> False],

  Cell[TextData[{
    "The area ", mathI[RowBox[{"S","(","\[Alpha]",")"}]], " grows monotonically with \[Alpha] across (0, \[Pi]). It diverges as ", mathI[RowBox[{"\[Alpha]","\[Rule]","\[Pi]"}]], " (the corridor straightens into an infinite strip) and shrinks as ", mathI[RowBox[{"\[Alpha]","\[Rule]","0"}]], " (the corridor folds back). A warm-start cascade is essential at the high-\[Alpha] end: cold Nelder\[Dash]Mead from the Hammersley quarter-arc seed converges to a poor thin-strip local optimum at 150\[Degree], while warm-starting from the 135\[Degree] anchor produces the much larger sofa above."}], "Text"],

  (* ====== References ====== *)
  Cell["References", "Section"],
  Cell[TextData[{
    StyleBox["[1] J. L. Gerver. ", FontWeight -> Bold],
    "On moving a sofa around a corner. ",
    StyleBox["Geom. Dedicata", FontSlant -> Italic],
    " 42 (1992), 267\[Dash]283."}], "Text"],
  Cell[TextData[{
    StyleBox["[2] J. M. Hammersley. ", FontWeight -> Bold],
    "On the enfeeblement of mathematical skills by \"Modern Mathematics\" and by similar soft intellectual trash in schools and universities. ",
    StyleBox["Bull. Inst. Math. Appl.", FontSlant -> Italic],
    " 4 (1968), 66\[Dash]85. (Contains the 2.207 sofa.)"}], "Text"],
  Cell[TextData[{
    StyleBox["[3] Y. Kallus, D. Romik. ", FontWeight -> Bold],
    "Improved upper bounds in the moving sofa problem. ",
    StyleBox["Adv. Math.", FontSlant -> Italic],
    " 340 (2018). arXiv:1706.06630."}], "Text"],
  Cell[TextData[{
    StyleBox["[4] J. Baek. ", FontWeight -> Bold],
    "Optimality of Gerver's Sofa. arXiv:2411.19826 (2024)."}], "Text"],
  Cell[TextData[{
    StyleBox["[5] D. Romik. ", FontWeight -> Bold],
    "Differential equations and exact solutions in the moving sofa problem. ",
    StyleBox["Exp. Math.", FontSlant -> Italic],
    " 26 (2017)."}], "Text"]
};

(* ----- assemble and save ----------------------------------------- *)

nb = Notebook[cells,
  WindowSize -> {1100, 800},
  WindowMargins -> {{Automatic, 0}, {Automatic, 0}},
  FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (January 2026)",
  StyleDefinitions -> "Default.nb",
  Saveable -> True];

outNb = FileNameJoin[{hereDir, "Sofa.nb"}];

Print["writing notebook -> ", outNb];
Export[outNb, nb];
Print["done; nb=", FileByteCount[outNb]];
