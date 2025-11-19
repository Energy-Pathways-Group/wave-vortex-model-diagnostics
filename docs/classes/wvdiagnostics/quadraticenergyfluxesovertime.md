---
layout: default
title: quadraticEnergyFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 92
mathjax: true
---

#  quadraticEnergyFluxesOverTime

Compute forcing fluxes over time.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxesOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) indices for time selection (default: Inf)

## Returns
+ `forcing_fluxes`  struct array with fluxes over time
+ `t`  Summary table of energy flux diagnostics

## Discussion

  Compute forcing fluxes over time
  Returns the energy fluxes from external forcing for each reservoir as a function of time.
 
              
