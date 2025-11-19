---
layout: default
title: plotEnergyTriadFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 76
mathjax: true
---

#  plotEnergyTriadFluxOverTime

Plot inertial flux for each reservoir over time


---

## Declaration
```matlab
 fig = plotInertialFluxOverTime(self,options)
```
## Parameters
+ `energyReservoirs`  vector of EnergyReservoir objects (default: [geostrophic, wave, total])
+ `visible`  figure visibility (default: "on")
+ `filter`  function handle to filter fluxes (default: @(v) v)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plots the energy flux between reservoirs due to inertial interactions as a function of time.
 
            
