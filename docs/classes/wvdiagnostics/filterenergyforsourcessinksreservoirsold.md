---
layout: default
title: filterEnergyForSourcesSinksReservoirsOld
parent: WVDiagnostics
grand_parent: Classes
nav_order: 46
mathjax: true
---

#  filterEnergyForSourcesSinksReservoirsOld

This function returns values assuming three reservoirs: geo, wave, and.


---

## Declaration
```matlab
 [sources, sinks, inertial, ddt, energy] = filterEnergyForSourcesSinksReservoirsOld(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `customNames`  (optional) input argument `customNames` (default: configureDictionary("string","string"))
+ `fluxTolerance`  (optional) input argument `fluxTolerance` (default: 1e-2)
+ `timeIndices`  (optional) indices specifying which time indices to use (default: Inf)
+ `shouldShowReservoirEnergy`  (optional) energy or enstrophy reservoir selection (default: true)
+ `shouldShowExactValues`  (optional) input argument `shouldShowExactValues` (default: true)
+ `shouldSeparateClosureRegion`  (optional) input argument `shouldSeparateClosureRegion` (default: true)

## Returns
+ `sources`  output value `sources`
+ `sinks`  output value `sinks`
+ `inertial`  output value `inertial`
+ `ddt`  output value `ddt`
+ `energy`  diagnosed energy as a function of time and/or scale

## Discussion

  This function returns values assuming three reservoirs: geo, wave, and
  damping. The damping resevoir is just scales below a threshold, wave or
  geostrophic. It also returns the exact and exact-damp resevoirs.
  The forcing struct has the the forcing on the three/two different
  reservoirs
  The inertial struct has the flux from the two reservoirs (wave
  geostrophic) to each other and to the damping region.
  The forcing struct also include the nonlinear advection, which has the
  flux to the damping region.
  The ddt struct contains the change in total energy, closing the energy
  budget.
 
                            
