---
layout: default
title: plotEnstrophyTriadFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 85
mathjax: true
---

#  plotEnstrophyTriadFluxOverTime

Plot enstrophy triad fluxes over time.


---

## Declaration
```matlab
 fig = plotEnstrophyTriadFluxOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `triadComponents`  (optional) vector of TriadFlowComponent objects to include (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
+ `timeIndices`  (optional) Indices of model times to use (default: Inf -> all times)
+ `visible`  (optional) Figure visibility (default: "on")
+ `filter`  (optional) Function handle accepting (v,t) used to preprocess each flux series before plotting (default: @(v,t) v)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot enstrophy triad fluxes over time.
  Plots time series of enstrophy triad fluxes computed with the quadratic
  approximation for the specified triad components. An optional filter may
  be applied to each time series before plotting. Axis labels use the class
  scaling properties (tscale, tscale_units, z_flux_scale, z_flux_scale_units).
 
                
