---
layout: default
title: plotSourcesSinksForReservoirGroup
parent: WVDiagnostics
grand_parent: Classes
nav_order: 92
mathjax: true
---

#  plotSourcesSinksForReservoirGroup

Plot sources, sinks, and reservoirs diagram for a reservoir group.


---

## Declaration
```matlab
 [fig, boxDiagram] = plotSourcesSinksForReservoirGroup(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `name`  (optional) Name of reservoir group to plot (default: "reservoir-damped-wave-geo")
+ `customNames`  (optional) dictionary mapping internal names to display names (default: configureDictionary("string","cell"))
+ `customColors`  (optional) dictionary mapping roles/names to RGB colors (default: dictionary with "source" and "sink" entries) (default: dictionary(["source","sink"], {[191 191 250]/255,[245 194 193]/255}))
+ `customReservoirOrder`  cell/array specifying reservoir ordering for layout (default: use group order)
+ `customForcing`  list of forcing names to include (default: all)
+ `fluxTolerance`  (optional) Minimum flux magnitude (in class flux units) to display arrows/labels (default: 1e-2)
+ `shouldShowUnits`  (optional) (optional, logical) Include units in labels (default: true)
+ `timeIndices`  (optional) Time indices to average over (default: Inf -> all times)
+ `shouldShowReservoirEnergy`  (optional) (optional, logical) Show reservoir energy as sublabels (default: true)
+ `shouldShowExactValues`  (optional) (optional, logical) Show exact total forcing values as sublabels (default: true)
+ `title`  (optional) Diagram title (default: "Energy Pathways")
+ `visible`  (optional) Figure visibility (default: "on")

## Returns
+ `fig`  handle to the generated figure
+ `boxDiagram`  BoxDiagram object used to draw the diagram (can be re-used or modified)

## Discussion

  Plot sources, sinks, and reservoirs diagram for a reservoir group.
  Generate a box-and-arrow diagram showing energy sources, sinks, reservoirs,
  and fluxes between them for a named reservoir group stored in the diagnostics
  NetCDF. Reads reservoir, forcing and inertial flux summaries for the requested
  timeIndices, formats labels and units using class scaling properties, and
  returns a drawn figure and the underlying BoxDiagram object.
 
                                  
