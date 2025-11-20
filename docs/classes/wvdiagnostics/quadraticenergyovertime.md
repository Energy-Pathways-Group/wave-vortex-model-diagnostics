---
layout: default
title: quadraticEnergyOverTime
parent: WVDiagnostics
grand_parent: Classes
nav_order: 100
mathjax: true
---

#  quadraticEnergyOverTime

Compute energy for each reservoir over time.


---

## Declaration
```matlab
 [reservoirs, t] = quadraticEnergyOverTime(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
+ `shouldIncludeExactTotalEnergy`  (optional) include exact total energy (default: true) (default: false)
+ `timeIndices`  (optional) indices for time selection (default: Inf)

## Returns
+ `reservoirs`  struct array with energy for each reservoir
+ `t`  time vector

## Discussion

  Compute energy for each reservoir over time
  Returns the energy in each specified reservoir as a function of time.
 
                
