---
layout: default
title: plotSourcesSinksReservoirsDiagramWithClosureRegion
parent: WVDiagnostics
grand_parent: Classes
nav_order: 94
mathjax: true
---

#  plotSourcesSinksReservoirsDiagramWithClosureRegion

Plot sources, sinks, and reservoirs diagram.


---

## Declaration
```matlab
 fig = plotSourcesSinksReservoirsDiagram(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `energyReservoirs`  (optional) vector of EnergyReservoir objects (default: [geostrophic, wave]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave])
+ `customNames`  (optional) dictionary for custom names (default: configureDictionary("string","string"))
+ `fluxTolerance`  (optional) tolerance for displaying fluxes (default: 1e-2)
+ `shouldShowUnits`  (optional) show units in labels (default: true)
+ `timeIndices`  (optional) indices for time averaging (default: Inf)
+ `shouldShowReservoirEnergy`  (optional) show reservoir energy (default: true)
+ `shouldShowExactValues`  (optional) input argument `shouldShowExactValues` (default: true)
+ `shouldSeparateClosureRegion`  (optional) input argument `shouldSeparateClosureRegion` (default: true)
+ `title`  (optional) diagram title (default: "Energy Pathways")
+ `visible`  (optional) figure visibility (default: "on")

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot sources, sinks, and reservoirs diagram
  Generates a diagram showing energy sources, sinks, and reservoirs, including fluxes between them.
 
                            
