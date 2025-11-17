---
layout: default
title: plotEnergyFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 73
mathjax: true
---

#  plotEnergyFluxOverTime

Plot energy fluxes for reservoirs as a function of time.


---

## Declaration
```matlab
 fig = plotEnergyFluxOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `approximation`  (optional) {'exact','quadratic'} which approximation to use (default: 'exact')
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects to include (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) Indices of model times to plot/average (default: Inf -> all times)
+ `visible`  (optional) Figure visibility (default: "on")
+ `filter`  (optional) Function handle to apply to flux series before plotting (default: @(v) v)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot energy fluxes for reservoirs as a function of time.
  Draws time series of energy fluxes (exact or quadratic approximation) into
  specified energy reservoirs. Supports selecting reservoirs, time indices,
  applying a filter to series prior to plotting, and configuring figure visibility.
 
                  
