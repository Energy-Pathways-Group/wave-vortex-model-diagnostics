---
layout: default
title: sparsePseudoRadialAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 124
mathjax: true
---

#  sparsePseudoRadialAxis

Create Sparse Pseudo-Radial Wavenumber Axis and Bin Assignments.


---

## Declaration
```matlab
 [kp,bins_0,bins_pm] = sparsePseudoRadialAxis(self)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `kp`  output value `kp`
+ `bins_0`  output value `bins_0`
+ `bins_pm`  output value `bins_pm`

## Discussion

  This function is used by create1DMirrorFluxes to create an efficient pseudo-radial wavenumber axis.
 
  The bin assignments, e.g. `bins_0` and `bins_pm`, can be used to create sparse matrices for efficient binning operations.
 
  ```matlab
  valid = ~isnan(bins_0);
  S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
  F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
  ```
 
            
