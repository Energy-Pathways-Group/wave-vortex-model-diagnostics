---
layout: default
title: quadraticEnstrophyTriadFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 120
mathjax: true
---

#  quadraticEnstrophyTriadFluxes

Return the enstrophy flux from the forcing terms.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxes(options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `triadComponents`  (optional) input argument `triadComponents` (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])

## Returns
+ `inertial_fluxes`  diagnosed flux values

## Discussion

  Return the enstrophy flux from the forcing terms
  Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
 
          
