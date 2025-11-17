---
layout: default
title: DST2
parent: WVDiagnostics
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  DST2

SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix.


---

## Declaration
```matlab
 matrix = DST2(N)
```
## Parameters
+ `N`  input argument `N`

## Returns
+ `matrix`  output value `matrix`

## Discussion

  SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix
  Forward scaling: 2/N (no endpoint tweaks).
  With this choice, the inverse (DST-III) is the plain sine matrix with
  the HIGHEST-frequency (last) COLUMN halved.
  If X = SineTransformForwardMatrix_DST2(N) * x, then
  x = InverseSineTransformMatrix_DST2(N) * X.
  See also: InverseSineTransformMatrix_DST2
 
        
