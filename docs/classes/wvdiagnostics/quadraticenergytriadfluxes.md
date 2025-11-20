---
layout: default
title: quadraticEnergyTriadFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 105
mathjax: true
---

#  quadraticEnergyTriadFluxes

Return the energy flux from the inertial terms, specified as triad components.


---

## Declaration
```matlab
 inertial_fluxes = quadraticEnergyTriadFluxes(options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
+ `triadComponents`  (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]. (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])

## Returns
+ `inertial_fluxes`  an array of structs

## Discussion

  Return the energy flux from the inertial terms, specified as triad components
  Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
 
            
