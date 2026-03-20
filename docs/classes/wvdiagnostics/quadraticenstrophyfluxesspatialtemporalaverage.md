---
layout: default
title: quadraticEnstrophyFluxesSpatialTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 113
mathjax: true
---

#  quadraticEnstrophyFluxesSpatialTemporalAverage

Compute spatial-temporal average of the qgpv enstrophy fluxes.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnstrophyFluxesSpatialTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `enstrophy_fluxes`  diagnosed flux values

## Discussion

  Compute spatial-temporal average of the qgpv enstrophy fluxes
  Returns the spatial-temporal average of the qgpv enstrophy fluxes from external forcing
 
          
