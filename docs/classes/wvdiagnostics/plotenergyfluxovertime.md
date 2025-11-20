---
layout: default
title: plotEnergyFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 72
mathjax: true
---

#  plotEnergyFluxOverTime

Plot energy fluxes as a function of time.


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

  You can plot either the exact or quadratic approximations to the energy fluxes.
 
  If you plot the quadratic fluxes, you can specify which energy reservoirs to include, which will create a subplot for each reservoir.
 
                  
