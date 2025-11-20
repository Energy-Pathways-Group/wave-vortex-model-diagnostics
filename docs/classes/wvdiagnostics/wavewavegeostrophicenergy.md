---
layout: default
title: waveWaveGeostrophicEnergy
parent: WVDiagnostics
grand_parent: Classes
nav_order: 141
mathjax: true
---

#  waveWaveGeostrophicEnergy

Note that.


---

## Declaration
```matlab
 E0 = waveWaveGeostrophicEnergy(wvt,mask)
```
## Parameters
+ `wvt`  WVDiagnostics object
+ `mask`  input argument `mask`

## Returns
+ `E0`  output value `E0`

## Discussion

  Note that
  energy = waveWaveGeostrophicEnergy(wvt,1,1);
  should produce the same answer as
  flow = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
  [Fp_w,Fm_w,F0_w] = wvt.nonlinearFluxForFlowComponents(flow,flow);
  [Ep_w,Em_w,E0_w] = wvt.energyFluxFromNonlinearFlux(Fp_w,Fm_w,F0_w);
  sum(E0_w(:))
 
          
