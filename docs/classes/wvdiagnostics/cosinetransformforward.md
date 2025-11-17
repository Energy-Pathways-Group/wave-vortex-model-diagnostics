---
layout: default
title: CosineTransformForward
parent: WVDiagnostics
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  CosineTransformForward

Fast Discrete Cosine Transform (DCT-I).


---

## Declaration
```matlab
 [xbar, f] = CosineTransformForward( t, x, varargin )
```
## Parameters
+ `t`  input argument `t`
+ `x`  input argument `x`
+ `varargin`  input argument `varargin`

## Returns
+ `xbar`  output value `xbar`
+ `f`  output value `f`

## Discussion

  CosineTransformForward  Fast Discrete Cosine Transform (DCT-I)
  xbar is returned in the same units as x. This is the finite length
  definition of a Fourier transform.
  f is returned in units of cycles.
  The cosine series would have the following sum,
  x(t) = xbar(1)/2 + sum( xbar(i)*cos(i*pi/T) )
  So note that the first coefficient double the average of the function.
  From this, Parseval's theorem is,
  x_sum = (1/T)*(sum(x(2:end-1).*x(2:end-1))*dt+x(1)*x(1)*dt/2 + x(end)*x(end)*dt/2);
  S_sum = ( S(1)/2 + sum(S(2:end-1)) + 2*S(end))*df;
  Where we've taken care to integrate only to the endpoints (and not
  beyond) in the x_sum. The S_sum includes a correction for the Nyquist.
 
              
