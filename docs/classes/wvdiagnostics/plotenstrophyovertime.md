---
layout: default
title: plotEnstrophyOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 79
mathjax: true
---

#  plotEnstrophyOverTime

Plot enstrophy over time


---

## Declaration
```matlab
 fig = plotEnstrophyOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `visible`  (optional) Figure visibility (default: "on")
+ `timeIndices`  (optional) Indices of model times to plot/average (default: Inf -> all times)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plots the quadratic and APV enstrophy as a function of time.
 
            
