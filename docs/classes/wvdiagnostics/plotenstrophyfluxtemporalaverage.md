---
layout: default
title: plotEnstrophyFluxTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 79
mathjax: true
---

#  plotEnstrophyFluxTemporalAverage

Plot temporally averaged enstrophy flux diagnostics.


---

## Declaration
```matlab
 fig = plotEnstrophyFluxTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) Time indices to average over (default: Inf -> all times)
+ `approximation`  (optional) {'quadratic','exact'} approximation to use (default: 'exact')
+ `axes`  (optional) Plot projection; one of {'jk','j','jWavenumber','k','k-pseudo-isotropic'} (default: 'jk')
+ `filter`  (optional) Function handle applied to plotted data (default: @(v) v)
+ `colormap`  (optional) Colormap for 'jk' axes (default: WVDiagnostics.crameri('-bam'))
+ `visible`  (optional) Figure visibility (default: "on")
+ `overSaturationFactor`  (optional) Scalar or two-element colormap limits (default: 10)
+ `simpleName`  (optional) Cell array of simple name strings for flux panels (default: [])

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Compute and plot temporally averaged enstrophy flux diagnostics (exact or
  quadratic approximation) on a chosen projection. Supports 'jk' (j vs
  radial wavelength) pcolor plots, 1D summaries along j or k axes, and a
  pseudo-isotropic radial projection. Allows applying a filter to the data,
  customizing the colormap and saturation, and controlling figure visibility.
 
                        
