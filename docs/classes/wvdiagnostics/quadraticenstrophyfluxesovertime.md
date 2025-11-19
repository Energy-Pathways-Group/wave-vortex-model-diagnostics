---
layout: default
title: quadraticEnstrophyFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 109
mathjax: true
---

#  quadraticEnstrophyFluxesOverTime

Compute enstrophy fluxes over time.


---

## Declaration
```matlab
 enstrophy_fluxes = quadraticEnstrophyFluxesOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `enstrophy_fluxes`  struct array with averaged fluxes
+ `t`  Summary table of enstrophy flux diagnostics

## Discussion

  Compute enstrophy fluxes over time
  Returns the enstrophy fluxes from external forcing
 
            
