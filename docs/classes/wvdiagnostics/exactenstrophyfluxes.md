---
layout: default
title: exactEnstrophyFluxes
parent: WVDiagnostics
grand_parent: Classes
nav_order: 40
mathjax: true
---

#  exactEnstrophyFluxes

Return the available potential enstrophy flux from the forcing terms, [j kRadial t].


---

## Declaration
```matlab
 forcing_fluxes = exactEnstrophyFluxes(options)
```
## Parameters
+ `self`  WVDiagnostics object

## Returns
+ `enstrophy_fluxes`  an array of structs

## Discussion

  Return the available potential enstrophy flux from the forcing terms, [j kRadial t]
  Reads from the diagnostics file and returns an array of structs with
  fields name, fancyName, and a field for each energy reservoir with size
  [j kRadial t]. This includes the nonlinear advection term.
 
        
