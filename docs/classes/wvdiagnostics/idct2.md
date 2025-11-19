---
layout: default
title: iDCT2
parent: WVDiagnostics
grand_parent: Classes
nav_order: 55
mathjax: true
---

#  iDCT2

InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix.


---

## Declaration
```matlab
 matrix = iDCT2(N)
```
## Parameters
+ `N`  input argument `N`

## Returns
+ `matrix`  output value `matrix`

## Discussion

  InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix
  This inverts CosineTransformForwardMatrix_DCT2 when the forward uses the
  2/N scaling above (with no extra row halving).
  Implementation: DCT-III with DC (first) COLUMN halved.
 
        
