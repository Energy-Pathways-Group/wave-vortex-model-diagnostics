---
layout: default
title: plotEnergyOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 76
mathjax: true
---

#  plotEnergyOverTime

Plot energy for each reservoir over time.


---

## Declaration
```matlab
 fig = plotEnergyOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
+ `shouldIncludeExactTotalEnergy`  (optional) include exact total energy (default: true)
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)
+ `visible`  (optional) figure visibility (default: "on")

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot energy for each reservoir over time
  Plots the energy in each specified reservoir as a function of time.
 
                
