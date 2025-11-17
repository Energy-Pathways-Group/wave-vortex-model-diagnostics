---
layout: default
title: createReservoirGroup
parent: WVDiagnostics
grand_parent: Classes
nav_order: 27
mathjax: true
---

#  createReservoirGroup

Create a reservoir group in the diagnostics NetCDF and populate it.


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

  Create a reservoir group in the diagnostics NetCDF and populate it.
  Create (or open) a diagnostics group containing reservoir definitions as defined by specified flow components.
  The function will create variables and dimensions as needed, then loop over
  the requested time indices to compute and write reservoir diagnostics.
 
              
