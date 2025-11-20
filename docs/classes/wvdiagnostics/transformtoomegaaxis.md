---
layout: default
title: transformToOmegaAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 133
mathjax: true
---

#  transformToOmegaAxis

transforms in the from (j,kRadial) to omegaAxis.


---

## Declaration
```matlab
 [varargout] = transformToRadialWavenumber(varargin)
```
## Parameters
+ `self`  WVDiagnostics object
+ `varargin`  variables with dimensions $$(j,kl)$$

## Returns
+ `varargout`  variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

## Discussion

  transforms in the from (j,kRadial) to omegaAxis
  Sums all the variance/energy in radial bins `kPseudoRadial`.
 
          
