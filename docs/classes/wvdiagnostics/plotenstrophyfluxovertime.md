---
layout: default
title: plotEnstrophyFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 77
mathjax: true
---

#  plotEnstrophyFluxOverTime

Plot enstrophy fluxes for reservoirs as a function of time.


---

## Declaration
```matlab
 fig = plotEnstrophyFluxOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `approximation`  (optional) {'quadratic','exact'} which approximation to use (default: 'exact')
+ `timeIndices`  (optional) Indices of model times to plot/average (default: Inf -> all times)
+ `visible`  (optional) Figure visibility (default: "on")
+ `filter`  (optional) Function handle to filter flux time series before plotting (default: @(v) v)
+ `shouldShowNonlinearAdvection`  (optional, logical) Show nonlinear advection term (default: true)
+ `shouldShowTotal`  (optional, logical) Show summed total flux (default: true)
+ `shouldShowDtEnstrophy`  (optional, logical) Show dZ/dt series for comparison (default: true)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Draws time series of enstrophy fluxes (exact or quadratic approximation)
  into specified enstrophy reservoirs. Supports selecting approximation,
  time indices, simple filtering for visualisation, and toggles for showing
  nonlinear advection, total flux, and dZ/dt.
 
  Some useful filters:
  filter=@(v,t) movmean(v,21);
  filter=@(v,t) cumtrapz(t,v)./(t+1)
 
                      
