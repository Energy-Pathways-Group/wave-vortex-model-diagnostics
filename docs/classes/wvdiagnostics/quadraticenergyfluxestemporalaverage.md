---
layout: default
title: quadraticEnergyFluxesTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 94
mathjax: true
---

#  quadraticEnergyFluxesTemporalAverage

Compute temporally averaged forcing fluxes.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxesTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `forcing_fluxes`  struct array with averaged fluxes

## Discussion

  Compute temporally averaged forcing fluxes
  Returns the temporally averaged energy fluxes from external forcing for each reservoir.
 
            
