---
layout: default
title: quadraticEnergyPrimaryTriadFluxesTemporalAverage2D
parent: WVDiagnostics
grand_parent: Classes
nav_order: 109
mathjax: true
---

#  quadraticEnergyPrimaryTriadFluxesTemporalAverage2D

outputGrid determines whether or not the fluxes get downsampled to the.


---

## Declaration
```matlab
 [inertial_fluxes_g, inertial_fluxes_w, k, j] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)
+ `outputGrid`  (optional) input argument `outputGrid` (default: "full")

## Returns
+ `inertial_fluxes_g`  diagnosed flux values
+ `inertial_fluxes_w`  diagnosed flux values
+ `k`  output value `k`
+ `j`  output value `j`

## Discussion

  outputGrid determines whether or not the fluxes get downsampled to the
  sparse grid, or up-sampled to the full grid.
 
                  
