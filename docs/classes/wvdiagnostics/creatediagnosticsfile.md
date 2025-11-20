---
layout: default
title: createDiagnosticsFile
parent: WVDiagnostics
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  createDiagnosticsFile

Create a new diagnostics file and compute diagnostics from WVModel output.


---

## Declaration
```matlab
 createDiagnosticsFile(self,options)
```
## Parameters
+ `self`  WVDiagnostics object.
+ `stride`  (optional) Stride for time sampling (default: 1).
+ `timeIndices`  Indices of time steps to process.
+ `filename`  Output path for the diagnostics file.
+ `shouldMeasureAntialiasingFlux`  (optional) (optional, logical) If true, computes antialiasing flux diagnostics (default: false).
+ `shouldUseHigherOrderFlux`  (optional) (optional, logical) If true, uses higher-order flux calculation (default: false).

## Discussion

  Create a new diagnostics file and compute diagnostics from WVModel output.
  Initializes a diagnostics NetCDF file, computes and stores energy, enstrophy, APE, APV, and flux diagnostics for each time step.
 
                
