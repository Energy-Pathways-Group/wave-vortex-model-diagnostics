---
layout: default
title: quadraticEnergyFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 96
mathjax: true
---

#  quadraticEnergyFluxes

Return the energy flux from the forcing terms.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxes(options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])

## Returns
+ `forcing_fluxes`  an array of structs

## Discussion

  Return the energy flux from the forcing terms
  Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
 
          
