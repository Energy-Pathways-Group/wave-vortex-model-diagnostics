---
layout: default
title: sparseJKAxisBinMatrices
parent: WVDiagnostics
grand_parent: Classes
nav_order: 120
mathjax: true
---

#  sparseJKAxisBinMatrices

Bin matrices for sparseKRadialAxis and sparseJWavenumberAxis


---

## Declaration
```matlab
 [S_0, S_pm, mask_0, mask_pm] = sparseJKAxisBinMatrices(self)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `S_0`  output value `S_0`
+ `S_pm`  output value `S_pm`
+ `mask_0`  output value `mask_0`
+ `mask_pm`  output value `mask_pm`

## Discussion

  Used by create2DMirrorFluxes to create sparse binning matrices for the (j,k) axes.
 
              
