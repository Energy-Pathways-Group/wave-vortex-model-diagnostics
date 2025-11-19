---
layout: default
title: SineTransformForward
parent: WVDiagnostics
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  SineTransformForward

Fast Discrete Sine Transform (DST-I).


---

## Parameters
+ `t`  input argument `t`
+ `x`  input argument `x`
+ `varargin`  input argument `varargin`

## Returns
+ `xbar`  output value `xbar`
+ `f`  output value `f`

## Discussion

  SineTransformForward  Fast Discrete Sine Transform (DST-I)
  xbar is returned in the same units as x. This is the finite length
  definition of a Fourier transform.
  f is returned in units of cycles, and only contains the potentially
  nonzero frequency i.e., the zero and nyquist are not included, because
  they must be zero.
  The following relationship is satisfied:
  [f, xbar] = SineTransformForward(t,x);
  S = T*(xbar .* conj(xbar));
  x_sum = (1/T)*sum(x.*x)*dt;
  S_sum = sum(S)*df;
  x_sum == S_sum
  Using,
  N=33; % total points
  T=1.0; % total time length
  t=T*(0:(N-1))'/(N-1);
  This definition of t *includes* the end points which *must* be zero.
 
              
