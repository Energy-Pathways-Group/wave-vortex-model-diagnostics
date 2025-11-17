---
layout: default
title: plotEnergyTriadFluxOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 79
mathjax: true
---

#  plotEnergyTriadFluxOverTime

Plot inertial flux for each reservoir over time.


---

## Declaration
```matlab
 fig = plotInertialFluxOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave])
+ `triadComponents`  (optional) input argument `triadComponents` (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)
+ `visible`  (optional) figure visibility (default: "on")
+ `filter`  (optional) function handle to filter fluxes (default: @(v) v)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot inertial flux for each reservoir over time
  Plots the energy flux between reservoirs due to inertial interactions as a function of time.
 
                  
