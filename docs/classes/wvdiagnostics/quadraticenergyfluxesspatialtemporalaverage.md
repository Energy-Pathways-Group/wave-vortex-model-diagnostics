---
layout: default
title: quadraticEnergyFluxesSpatialTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 98
mathjax: true
---

#  quadraticEnergyFluxesSpatialTemporalAverage

Compute spatial-temporal average of forcing fluxes.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxesSpatialTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `forcing_fluxes`  struct array with averaged fluxes

## Discussion

  Compute spatial-temporal average of forcing fluxes
  Returns the spatial-temporal average of energy fluxes from external forcing for each reservoir.
 
            
