---
layout: default
title: exactEnergyFluxesTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 37
mathjax: true
---

#  exactEnergyFluxesTemporalAverage

Temporally averaged exact energy fluxes.


---

## Declaration
```matlab
 enstrophy_fluxes = exactEnergyFluxesTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `forcing_fluxes`  diagnosed flux values

## Discussion

  Temporally averaged exact energy fluxes
  Returns the temporally averaged enstrophy fluxes from external forcing for each reservoir.
 
          
