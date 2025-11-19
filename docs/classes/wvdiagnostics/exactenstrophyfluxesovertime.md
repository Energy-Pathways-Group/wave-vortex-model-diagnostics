---
layout: default
title: exactEnstrophyFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 40
mathjax: true
---

#  exactEnstrophyFluxesOverTime

Compute exact enstrophy fluxes over time.


---

## Declaration
```matlab
 forcing_fluxes = exactEnstrophyFluxesOverTime(self)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)

## Returns
+ `forcing_fluxes`  struct array with exact fluxes
+ `t`  Summary table of enstrophy flux diagnostics

## Discussion

  Compute exact enstrophy fluxes over time
  Returns the exact enstrophy fluxes from external forcing for each time step.
 
            
