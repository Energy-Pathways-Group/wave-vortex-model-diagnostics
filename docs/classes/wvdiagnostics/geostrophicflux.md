---
layout: default
title: geostrophicFlux
parent: WVDiagnostics
grand_parent: Classes
nav_order: 55
mathjax: true
---

#  geostrophicFlux

Compute the spatial and spectral fluxes for constant stratification


---

## Declaration
```matlab
 [spatialFlux, spectralFlux] = geostrophicFlux(self)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `spatialFlux`  struct containing `ggg`, `ggw`, `ggw_tx`, `wwg_tx`
+ `spectralFlux`  struct containing `ggg`, `ggw`, `ggw_tx`, `wwg_tx`

## Discussion

  This is from the source-vector form of the fluxes, derived by Jonathan
  Lilly and Jeffrey Early, and is only valid for constant stratification at
  the moment.
 
 
          
