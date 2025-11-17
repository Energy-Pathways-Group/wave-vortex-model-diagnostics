---
layout: default
title: quadraticEnstrophyFluxesTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 118
mathjax: true
---

#  quadraticEnstrophyFluxesTemporalAverage

Compute temporally averaged enstrophy fluxes.


---

## Declaration
```matlab
 enstrophy_fluxes = quadraticEnstrophyFluxesTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `enstrophy_fluxes`  struct array with averaged fluxes

## Discussion

  Compute temporally averaged enstrophy fluxes
  Returns the temporally averaged enstrophy fluxes from external forcing for each reservoir.
 
          
