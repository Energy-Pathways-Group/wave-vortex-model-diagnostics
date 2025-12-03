---
layout: default
title: addTriadFluxesForReservoirGroupAtTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  addTriadFluxesForReservoirGroupAtTime

Add Triad Fluxes For Reservoir Group At Time.


---

## Declaration
```matlab
 addTriadFluxesForReservoirGroupAtTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `triadVar`  input argument `triadVar`
+ `flowComponents`  input argument `flowComponents`
+ `wvt`  WVDiagnostics object
+ `outputIndex`  input argument `outputIndex`

## Discussion

  Compute triad energy-flux contributions for a reservoir group at the
  specified time index and write the integrated per-triad energy
  transfers into the corresponding `triadVar("T_i_j_k")` variables.
 
              
