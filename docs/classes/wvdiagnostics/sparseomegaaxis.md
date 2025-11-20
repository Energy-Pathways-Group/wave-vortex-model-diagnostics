---
layout: default
title: sparseOmegaAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 123
mathjax: true
---

#  sparseOmegaAxis

Sparse Omega Axis.


---

## Declaration
```matlab
 [omegaAxis,bins_omega] = sparseOmegaAxis(self)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `omegaAxis`  output value `omegaAxis`
+ `bins_omega`  output value `bins_omega`

## Discussion

  This function is used by create1DMirrorFluxes to create an efficient omega axis.
 
  The bin assignments, `bins_omega`, can be used to create sparse matrices for efficient binning operations.
 
  ```matlab
  valid = ~isnan(bins_0);
  S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
  F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
  ```
 
          
