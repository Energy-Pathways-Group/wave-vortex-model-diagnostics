---
layout: default
title: exactEnergyFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 36
mathjax: true
---

#  exactEnergyFluxesOverTime

Exact energy fluxes over time.


---

## Declaration
```matlab
 forcing_fluxes = exactEnergyFluxesOverTime(self)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)

## Returns
+ `energy_fluxes`  diagnosed flux values
+ `t`  Summary table of energy flux diagnostics

## Discussion

  Exact energy fluxes over time
  Returns the exact energy fluxes from external forcing for each time step.
 
            
