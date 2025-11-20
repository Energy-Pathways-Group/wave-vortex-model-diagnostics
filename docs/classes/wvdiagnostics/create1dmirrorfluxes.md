---
layout: default
title: create1DMirrorFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  create1DMirrorFluxes

Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.


---

## Declaration
```matlab
 create1DMirrorFluxes(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `stride`  (optional) Stride for time sampling (default: 1).
+ `timeIndices`  Indices of time steps to process (default: all times in diagnostics file).

## Discussion

  Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
  Adds 1D mirror-flux variables (kp, omegaAxis, kePeAxis) to an existing diagnostics
  NetCDF file and computes mirror flux summaries (e.g. F_wwg_kp, pi_w_wwg_kp,
  F_ggw_kp, pi_g_ggw_kp, and the corresponding omega/kePe projections) from the
  current WVTransform for the requested time indices. The function will create
  required dimensions/variables if they do not already exist and then populate
  the variables for each requested time index.
 
          
