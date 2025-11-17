---
layout: default
title: quadraticEnstrophyTriadFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 121
mathjax: true
---

#  quadraticEnstrophyTriadFluxesOverTime

Compute enstrophy inertial (aka, triad) fluxes over time.


---

## Declaration
```matlab
 enstrophy_fluxes = quadraticEnstrophyFluxesOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `triadComponents`  (optional) input argument `triadComponents` (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
+ `timeIndices`  (optional) indices for time averaging (default: Inf)

## Returns
+ `enstrophy_fluxes`  struct array with averaged fluxes
+ `t`  Summary table of enstrophy flux diagnostics

## Discussion

  Compute enstrophy inertial (aka, triad) fluxes over time
  Returns the enstrophy fluxes from external forcing
 
              
