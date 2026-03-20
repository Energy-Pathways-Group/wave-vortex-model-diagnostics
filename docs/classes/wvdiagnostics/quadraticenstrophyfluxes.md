---
layout: default
title: quadraticEnstrophyFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 111
mathjax: true
---

#  quadraticEnstrophyFluxes

Return the enstrophy flux from the forcing terms.


---

## Declaration
```matlab
 forcing_fluxes = quadraticEnergyFluxes(options)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `enstrophy_fluxes`  diagnosed flux values

## Discussion

  Return the enstrophy flux from the forcing terms
  Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
 
        
