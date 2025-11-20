---
layout: default
title: sparseKePeAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 123
mathjax: true
---

#  sparseKePeAxis

Sparse kenetic-potential energy ratio Axis.


---

## Declaration
```matlab
 [kePeAxis,bins_kepe] = sparseKePeAxis(self)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `kePeAxis`  output value `kePeAxis`
+ `bins_kepe`  output value `bins_kepe`

## Discussion

  This function is used by create1DMirrorFluxes to create an efficient ke-pe axis.
 
  The bin assignments, `bins_kepe`, can be used to create sparse matrices for efficient binning operations.
 
  ```matlab
  valid = ~isnan(bins_0);
  S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
  F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
  ```
 
          
