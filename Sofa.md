```wl
In[]:= reg = RegionUnion[DiscretizeGraphics[Graphics@Rectangle[{-5, -1}, {2, 1}]], DiscretizeGraphics[Graphics@Rectangle[{0, -1}, {2, -6}]], AspectRatio -> 1]
```

![0qp7p8hx2gll3](img/0qp7p8hx2gll3.png)

```wl
In[]:= Manipulate[RegionIntersection[reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{x, 0}]],RotationTransform[alpha, {0, -1}]], PlotRange -> {{-6, 6}, {-6, 6}}], {x, 0, 2}, {alpha, 0, 2 Pi}]
```

MarkdownTools`Private`ExportAnimatedImage[-AnimatedImage-, {GeneratedAssetLocation -> /Users/thiel/GitHub/SofaProblem/img/, WolframLanguageTag -> wl}]

![05wpg08sv6lnj](img/05wpg08sv6lnj.png)

```wl
In[]:= Manipulate[RegionPlot[{reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{x, 0}]], RotationTransform[-alpha, {0, -1}]]}, PlotRange -> {{-6, 6}, {-6, 6}}], {x, 0, 2}, {alpha, 0, 2 Pi}]
```

MarkdownTools`Private`ExportAnimatedImage[-AnimatedImage-, {GeneratedAssetLocation -> /Users/thiel/GitHub/SofaProblem/img/, WolframLanguageTag -> wl}]

![1fayrjftngf0d](img/1fayrjftngf0d.png)

```wl
In[]:= Manipulate[RegionPlot[{reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{-2 Cos[alpha] + 1, -1 + 2 Cos[alpha]}]], RotationTransform[-alpha + Pi/4, {0, -1}]]}, PlotRange -> {{-6, 6}, {-6, 6}}], {alpha, 0, Pi/2}]
```

MarkdownTools`Private`ExportAnimatedImage[-AnimatedImage-, {GeneratedAssetLocation -> /Users/thiel/GitHub/SofaProblem/img/, WolframLanguageTag -> wl}]

![1kz46mh3rz4cd](img/1kz46mh3rz4cd.png)

```wl
In[]:= \[AliasDelimiter]\[AliasDelimiter]
```

```wl
In[]:= TransformedRegion[reg, RotationTransform[2.]]
```

![0kkq2tm8lsewr](img/0kkq2tm8lsewr.png)

```wl
In[]:= RegionIntersection[reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{x, 0}]],RotationTransform[alpha]]]
```

![0aw89cygi1guq](img/0aw89cygi1guq.png)

```wl
In[]:= ListLinePlot[Table[{1. Cos[k alpha], 1. Sin[k alpha]} /. k -> 2., {alpha, -Pi/4, Pi/4, (Pi/2.)/100.}], AspectRatio -> 1]
```

```wl
In[]:= 
```

```wl
In[]:= regs = Flatten[Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{0.01 Cos[k alpha], 0.01 Sin[k alpha]}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2., (Pi/2.)/10.}, {k, 1, 9}]];
```

```wl
In[]:= regs
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1h242zm5ut4c3](img/1h242zm5ut4c3.png)

```wl
In[]:= \[AliasDelimiter]\[AliasDelimiter]
```

```wl
In[]:= RegionIntersection @@ Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{Cos[10 alpha], Sin[10 alpha]}]], RotationTransform[0.1*alpha]], {alpha, 0, Pi/2., 0.1}]
```

![0d772lcsxpg4d](img/0d772lcsxpg4d.png)

```wl
In[]:= (*Alpha needs to go from 0,Pi/2*)
```

```wl
In[]:= ListLinePlot[Accumulate@Table[{2/100. Cos[k alpha], 2/100. Sin[k alpha]} /. alpha -> 0.6, {k, 1, 100, 1}], AspectRatio -> 1]
```

![067aca0178o0b](img/067aca0178o0b.png)

```wl
In[]:= ListPlot[Table[{-1 Cos[alpha], 1}, {alpha, 0, Pi/2.}]]
```

```wl
In[]:= {-1, 1}, {1, -1}
```

```wl
In[]:= Plot[{-Cos[2 alpha], Cos[2 alpha]}, {alpha, 0, Pi/2}, AspectRatio -> 1]
```

![0obaoystjuaix](img/0obaoystjuaix.png)

```wl
In[]:= Plot[{-2 Cos[alpha] + 1, -1 + 2 Cos[alpha]}, {alpha, 0, Pi/2}, AspectRatio -> 1]
```

![0zq7jo4d8jp9i](img/0zq7jo4d8jp9i.png)

```wl
In[]:= Plot[{-Cos[2 alpha], Cos[2 alpha]}, {alpha, 0, Pi/2}, AspectRatio -> 1]
```

![1uuktlvmtziy7](img/1uuktlvmtziy7.png)

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-2 Cos[alpha] + 1, -1 + 2 Cos[alpha]}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.03}];
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-Cos[2 alpha], Cos[2 alpha]}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![0vmyzrsyp52ft](img/0vmyzrsyp52ft.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.18185
```

```wl
In[]:= %/4
```

```wl
Out[]= 2.04546
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-0.9 Cos[2 alpha], 0.9 Cos[2 alpha]}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![16mq9yb06lprz](img/16mq9yb06lprz.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.18533
```

```wl
In[]:= Plot[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.5, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.5}, {alpha, 0, Pi/2}]
```

![0rokjx6bv8jjc](img/0rokjx6bv8jjc.png)

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1ous2fkcy4iud](img/1ous2fkcy4iud.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.28222
```

```wl
In[]:= %/4
```

```wl
Out[]= 2.07055
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.04 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![0lw1ouxqy7j1l](img/0lw1ouxqy7j1l.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.26161
```

```wl
In[]:= regs // Length
```

```wl
Out[]= 79
```

```wl
In[]:= FoldList[RegionIntersection, regs[[1 ;; 79]]]
```

```wl
In[]:= reg = RegionUnion[DiscretizeGraphics[Graphics@Rectangle[{-7, -1}, {2, 1}]], DiscretizeGraphics[Graphics@Rectangle[{0, -1}, {2, -8}]], AspectRatio -> 1]
```

![19ddqx5h57jrc](img/19ddqx5h57jrc.png)

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1.1}]], RotationTransform[-alpha, {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= FoldList[RegionIntersection, regs]
```

```wl
In[]:= Manipulate[RegionPlot[{reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1}]], RotationTransform[-alpha + Pi/4., {0, -1}]]}, PlotRange -> {{-10, 10}, {-10, 10}}], {alpha, 0, Pi/2}]
```

MarkdownTools`Private`ExportAnimatedImage[-AnimatedImage-, {GeneratedAssetLocation -> /Users/thiel/GitHub/SofaProblem/img/, WolframLanguageTag -> wl}]

```wl
In[]:= Plot[{-alpha + Pi/4, Sign[-alpha + Pi/4.] Sqrt[Abs[(-alpha + Pi/4.)]*4/Pi]*Pi/4.}, {alpha, 0, Pi/2}]
```

![1iw0cx9koduks](img/1iw0cx9koduks.png)

```wl
In[]:= Manipulate[RegionPlot[{reg, TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1}]], RotationTransform[Sign[-alpha + Pi/4.] Sqrt[Abs[(-alpha + Pi/4.)]*4/Pi]*Pi/4., {0, -1}]]}, PlotRange -> {{-10, 10}, {-10, 10}}], {alpha, 0, Pi/2}]
```

MarkdownTools`Private`ExportAnimatedImage[-AnimatedImage-, {GeneratedAssetLocation -> /Users/thiel/GitHub/SofaProblem/img/, WolframLanguageTag -> wl}]

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1, 1 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.6*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![159vqk57v5j21](img/159vqk57v5j21.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.40544
```

```wl
In[]:= %/4
```

```wl
Out[]= 2.10136
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.05 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1, 1.05 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.6*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1a988d1cc0de9](img/1a988d1cc0de9.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.41839
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^1}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.66*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1n8v4ab73nqlm](img/1n8v4ab73nqlm.png)

```wl
In[]:= Area[%658]
```

```wl
Out[]= 8.45864
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.98, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.98}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.66*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![0hhtkf7gdiieq](img/0hhtkf7gdiieq.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.46654
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.94, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.94}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.66*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![07siry4q4a57z](img/07siry4q4a57z.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.47997
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.9, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.9}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.66*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1lodshi2oaxsy](img/1lodshi2oaxsy.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.4896
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.66*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.02}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![0yblmvnq6kcyd](img/0yblmvnq6kcyd.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.49304
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.63*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1lqew2kmqoea1](img/1lqew2kmqoea1.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.49835
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8, 1.01 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.8}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.6*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![15qdhd8xmlbdp](img/15qdhd8xmlbdp.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.50305
```

```wl
In[]:= %/4.
```

```wl
Out[]= 2.12576
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.2 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.7, 1.2 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.7}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.53*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1mrkoruesmrwy](img/1mrkoruesmrwy.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.55506
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68, 1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.53*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1l5hs1wwmqw2c](img/1l5hs1wwmqw2c.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.55843
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68, 1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.52*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![0eii1lqduc9yo](img/0eii1lqduc9yo.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.5583
```

```wl
In[]:= %/4.
```

```wl
Out[]= 2.13957
```

```wl
In[]:= regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68, 1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.52*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.01}];
```

```wl
In[]:= Fold[RegionIntersection, regs]
```

![1vbvnx49ktwnx](img/1vbvnx49ktwnx.png)

```wl
In[]:= Area[%]
```

```wl
Out[]= 8.51918
```

```wl
In[]:= %/4.
```

```wl
Out[]= 2.1298
```

```wl
In[]:= Plot[{-1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68, 1.225 Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^0.68}, {alpha, 0, Pi/2.}]
```

![0k2y23k67t3sn](img/0k2y23k67t3sn.png)

```wl
In[]:= \[AliasDelimiter]
```

```wl
In[]:= Plot[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^0.52*Pi/4., {alpha, 0, Pi/2.}]
```

![0itj2svnfq2ci](img/0itj2svnfq2ci.png)

```wl
In[]:= NMaximize
```

```wl
In[]:= results = {}; Monitor[Table[regs = Table[TransformedRegion[TransformedRegion[reg, TranslationTransform[{-amp Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^exp1, amp Sign[Cos[2 alpha]] Abs[Cos[2 alpha]]^exp1}]], RotationTransform[Sign[-alpha + Pi/4.] (Abs[(-alpha + Pi/4.)]*4/Pi)^exp2*Pi/4., {0, -1}]], {alpha, 0, Pi/2, 0.015}]; AppendTo[results, {amp, exp1, exp2, Area[Fold[RegionIntersection, regs]]}];Export["~/Desktop/results.mx", results], {amp, 1.1, 1.6, 0.1}, {exp1, 0.55, 0.7, 0.03}, {exp2, 0.45, 0.6, 0.03}];, {amp, exp1, exp2}]
```

![1cqguv7ke92du](img/1cqguv7ke92du.png)