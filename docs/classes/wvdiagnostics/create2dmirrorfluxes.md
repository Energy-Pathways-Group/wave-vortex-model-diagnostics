---
layout: default
title: create2DMirrorFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  create2DMirrorFluxes

Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.


---

## Declaration
```matlab
 create2DMirrorFluxes(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `stride`  (optional) Stride between model time steps when creating new output (default: 1)
+ `timeIndices`  Explicit list of model time indices to process (default: uses stride over all times)
+ `mirrorTriad`  (optional) Which mirror triad to compute, one of "wwg" or "ggw" (default: "wwg")
+ `shouldOverwriteExisting`  (optional) (optional, logical) If true, overwrite any existing group data rather than appending (default: false)

## Discussion

  Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
  Generates 2D mirror-flux summaries (triad primary and mirror variables)
  binned on the js/ks grid for the specified triad type ("wwg" or "ggw").
  Creates a diagnostics group (mirror-flux-2d-<triad>) and variables if they
  do not exist, then computes and appends the per-time-index 2D flux fields.
 
              
