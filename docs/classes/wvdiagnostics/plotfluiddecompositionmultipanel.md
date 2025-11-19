---
layout: default
title: plotFluidDecompositionMultipanel
parent: WVDiagnostics
grand_parent: Classes
nav_order: 82
mathjax: true
---

#  plotFluidDecompositionMultipanel

Plot a multipanel decomposition of the fluid state (wave vs geostrophic).


---

## Declaration
```matlab
 fig = plotFluidDecompositionMultipanel(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `visible`  (optional) Figure visibility (default: "on")
+ `iTime`  Time index to display; when provided sets self.iTime (default: self.iTime)
+ `title`  Figure title; default uses the model time in days
+ `yForXZSlice`  y coordinate (meters) at which to take the x-z slice; default uses center y

## Returns
+ `fig`  Handle to the generated figure

## Discussion

  Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
  Creates a compact multipanel figure showing horizontal (x-y) maps and
  vertical (x-z) sections of vertical vorticity for the total flow, the
  wave component, and the geostrophic component at a specified model time.
  Uses fields from the current WVTransform (self.wvt) and applies consistent
  colormap, plotting scales, and annotations. The x/y axes are presented in
  kilometers and the z (depth) axis in kilometers.
 
                
