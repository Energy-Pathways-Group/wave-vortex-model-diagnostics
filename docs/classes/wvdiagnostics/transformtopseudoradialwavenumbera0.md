---
layout: default
title: transformToPseudoRadialWavenumberA0
parent: WVDiagnostics
grand_parent: Classes
nav_order: 134
mathjax: true
---

#  transformToPseudoRadialWavenumberA0

transforms in the from (j,kRadial) to kPseudoRadial.


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

  transforms in the from (j,kRadial) to kPseudoRadial
  Sums all the variance/energy in radial bins `kPseudoRadial`.
 
          
