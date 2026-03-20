---
layout: default
title: waveWaveGeostrophicEnergyForMode
parent: WVDiagnostics
grand_parent: Classes
nav_order: 144
mathjax: true
---

#  waveWaveGeostrophicEnergyForMode

Note that.


---

## Declaration
```matlab
 E0 = waveWaveGeostrophicEnergyForMode(wvt,maskKU,maskKUx,Nj)
```
## Parameters
+ `wvt`  WVDiagnostics object
+ `maskKU`  input argument `maskKU`
+ `maskKUx`  input argument `maskKUx`
+ `Nj`  input argument `Nj`

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
 
              
