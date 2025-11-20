---
layout: default
title: exactEnergyFluxesSpatialTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 37
mathjax: true
---

#  exactEnergyFluxesSpatialTemporalAverage

Compute spatial-temporal average of the exact forcing fluxes.


---

## Declaration
```matlab
 forcing_fluxes = exactEnergyFluxesSpatialTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `forcing_fluxes`  struct array with averaged fluxes

## Discussion

  Compute spatial-temporal average of the exact forcing fluxes
  Returns the spatial-temporal average of the exact energy fluxes from external forcing
 
          
