---
layout: default
title: quadraticEnergyTriadFluxesTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 108
mathjax: true
---

#  quadraticEnergyTriadFluxesTemporalAverage

Computes the temporally averaged quadratic triad fluxes.


---

## Declaration
```matlab
 inertial_fluxes = quadraticEnergyTriadFluxesTemporalAverage(options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) indices specifying which time steps to average over. Defaults to Inf (all). (default: Inf)
+ `triadComponents`  (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]. (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])

## Returns
+ `inertial_fluxes`  an array of structs

## Discussion

  Computes the temporally averaged quadratic triad fluxes.
  Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial].
 
              
