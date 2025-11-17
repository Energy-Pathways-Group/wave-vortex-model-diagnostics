---
layout: default
title: geostrophicGeostrophicWaveEnergy
parent: WVDiagnostics
grand_parent: Classes
nav_order: 55
mathjax: true
---

#  geostrophicGeostrophicWaveEnergy

Note that.


---

## Declaration
```matlab
 Epm = geostrophicGeostrophicWaveEnergy(wvt,mask)
```
## Parameters
+ `wvt`  WVDiagnostics object
+ `mask`  input argument `mask`

## Returns
+ `Epm`  output value `Epm`

## Discussion

  Note that
  energy = waveWaveGeostrophicEnergy(wvt,1,1);
  should produce the same answer as
  flow = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
  [Fp_w,Fm_w,F0_w] = wvt.nonlinearFluxForFlowComponents(flow,flow);
  [Ep_w,Em_w,E0_w] = wvt.energyFluxFromNonlinearFlux(Fp_w,Fm_w,F0_w);
  sum(E0_w(:))
 
          
