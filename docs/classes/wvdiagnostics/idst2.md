---
layout: default
title: iDST2
parent: WVDiagnostics
grand_parent: Classes
nav_order: 57
mathjax: true
---

#  iDST2

InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix.


---

## Declaration
```matlab
 matrix = iDST2(N)
```
## Parameters
+ `N`  input argument `N`

## Returns
+ `matrix`  output value `matrix`

## Discussion

  InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix
  This inverts SineTransformForwardMatrix_DST2 when the forward uses the
  2/N scaling above (with no endpoint tweaks).
  Implementation: DST-III with the LAST column halved (k = N).
 
        
