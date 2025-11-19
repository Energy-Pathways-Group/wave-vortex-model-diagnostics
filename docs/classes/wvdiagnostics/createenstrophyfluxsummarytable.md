---
layout: default
title: createEnstrophyFluxSummaryTable
parent: WVDiagnostics
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  createEnstrophyFluxSummaryTable

Create a LaTeX table summarizing enstrophy source/sink fluxes.


---

## Declaration
```matlab
 tableString = createEnstrophyFluxSummaryTable(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `customNames`  (optional) dictionary mapping diagnostic names to custom labels (default: configureDictionary("string","string"))
+ `timeIndices`  (optional) time indices to average over (default: Inf -> all times)

## Returns
+ `tableString`  string containing the LaTeX table

## Discussion

  Create a LaTeX table summarizing enstrophy source/sink fluxes.
  Generates a LaTeX-formatted table listing enstrophy sources, sinks and totals
  using spatial-temporally averaged enstrophy flux diagnostics. The table
  reports both APV (exact) and quadratic (QGPV) flux measures and formats
  columns so ampersands are aligned for human-readable LaTeX output.
 
            
