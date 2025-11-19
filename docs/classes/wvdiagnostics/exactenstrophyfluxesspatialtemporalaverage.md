---
layout: default
title: exactEnstrophyFluxesSpatialTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 41
mathjax: true
---

#  exactEnstrophyFluxesSpatialTemporalAverage

Compute spatial-temporal average of the exact enstrophy fluxes.


---

## Declaration
```matlab
 forcing_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `enstrophy_fluxes`  diagnosed flux values

## Discussion

  Compute spatial-temporal average of the exact enstrophy fluxes
  Returns the spatial-temporal average of the exact enstrophy fluxes from external forcing
 
          
