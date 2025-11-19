---
layout: default
title: transformToPseudoRadialWavenumberApm
parent: WVDiagnostics
grand_parent: Classes
nav_order: 135
mathjax: true
---

#  transformToPseudoRadialWavenumberApm

transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial.


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

  transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial
  Sums all the variance/energy in radial bins `kPseudoRadial`.
 
          
