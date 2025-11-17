---
layout: default
title: DCT2
parent: WVDiagnostics
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  DCT2

CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix.


---

## Declaration
```matlab
 matrix = DCT2(N)
```
## Parameters
+ `N`  input argument `N`

## Returns
+ `matrix`  output value `matrix`

## Discussion

  CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix
  Forward scaling: 2/N (no extra endpoint halving here).
  With this choice, the inverse (DCT-III) is the plain cosine matrix with
  the DC (first) COLUMN halved.
  If X = CosineTransformForwardMatrix_DCT2(N) * x, then
  x = InverseCosineTransformMatrix_DCT2(N) * X.
  See also: InverseCosineTransformMatrix_DCT2
 
        
