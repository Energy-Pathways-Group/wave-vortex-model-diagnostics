---
layout: default
title: quadraticEnergyTriadFluxesOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 108
mathjax: true
---

#  quadraticEnergyTriadFluxesOverTime

Compute inertial fluxes over time.


---

## Declaration
```matlab
 inertial_fluxes = quadraticEnergyTriadFluxesOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
+ `triadComponents`  (optional) vector of TriadFlowComponent objects (default: [geostrophic_mda, wave]) (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)

## Returns
+ `inertial_fluxes`  struct array with fluxes over time
+ `t`  Summary table of energy flux diagnostics

## Discussion

  Compute inertial fluxes over time
  Returns the energy fluxes due to inertial interactions for each reservoir as a function of time.
 
                
