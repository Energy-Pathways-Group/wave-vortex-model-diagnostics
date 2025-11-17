---
layout: default
title: PoissonFlowFromFluxDCTI
parent: WVDiagnostics
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  PoissonFlowFromFluxDCTI

We will treat the first dimension as `x' and the second


---

## Discussion
dimension as `y'. This means that the flux in the usual form,
  which is j by kRadial, might need to be transposed to get
  what you want.
 
  [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
  quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
