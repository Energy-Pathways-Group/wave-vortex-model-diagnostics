---
layout: default
title: fluxesForReservoirGroup
parent: WVDiagnostics
grand_parent: Classes
nav_order: 49
mathjax: true
---

#  fluxesForReservoirGroup

Fluxes For Reservoir Group.


---

## Declaration
```matlab
 [transferFlux, forcingFlux, ddt, energy] = fluxesForReservoirGroup(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `outputfile`  file path or name
+ `name`  (optional) input argument `name` (default: "reservoir-damped-wave-geo")
+ `timeIndices`  indices specifying which time indices to use
+ `outputFormat`  (optional) input argument `outputFormat` (default: 'matrix')

## Returns
+ `transferFlux`  diagnosed flux values
+ `forcingFlux`  diagnosed flux values
+ `ddt`  output value `ddt`
+ `energy`  diagnosed energy as a function of time and/or scale

## Discussion

  fluxesForReservoirGroup is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
 
                      
