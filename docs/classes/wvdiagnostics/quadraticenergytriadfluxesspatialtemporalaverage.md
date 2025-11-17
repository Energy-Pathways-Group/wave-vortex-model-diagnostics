---
layout: default
title: quadraticEnergyTriadFluxesSpatialTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 113
mathjax: true
---

#  quadraticEnergyTriadFluxesSpatialTemporalAverage

Compute spatial-temporal average of inertial fluxes.


---

## Declaration
```matlab
 inertial_fluxes = quadraticEnergyTriadFluxesSpatialTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) indices for time averaging (default: Inf)
+ `triadComponents`  (optional) vector of TriadFlowComponent objects (default: [geostrophic_mda, wave]) (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])

## Returns
+ `inertial_fluxes`  struct array with averaged fluxes

## Discussion

  Compute spatial-temporal average of inertial fluxes
  Returns the spatial-temporal average of energy fluxes due to inertial interactions for each reservoir.
 
              
