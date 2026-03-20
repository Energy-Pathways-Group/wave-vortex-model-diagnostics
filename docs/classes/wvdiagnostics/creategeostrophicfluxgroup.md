---
layout: default
title: createGeostrophicFluxGroup
parent: WVDiagnostics
grand_parent: Classes
nav_order: 27
mathjax: true
---

#  createGeostrophicFluxGroup

Create a geostrophic flux group in the diagnostics NetCDF and populate it.


---

## Declaration
```matlab
 createReservoirGroup(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `outputfile`  NetCDFGroup to write into (default: self.diagfile).
+ `name`  (optional) Name of the reservoir group (default: "reservoir-damped-wave-geo").
+ `flowComponents`  WVFlowComponent array specifying reservoirs to create.
+ `timeIndices`  Time indices to process (default: all times in diagnostics file).

## Discussion

  Computes the geostrophic fluxes from geostrophicFlux, valid only for
  constant stratification..
 
              
