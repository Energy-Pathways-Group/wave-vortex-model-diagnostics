---
layout: default
title: plotEnergyFluxes1D
parent: WVDiagnostics
grand_parent: Classes
nav_order: 73
mathjax: true
---

#  plotEnergyFluxes1D

Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency.


---

## Declaration
```matlab
 fig = plotEnergyFluxes1D(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `timeIndices`  (optional) Indices of model times to average/plot (default: Inf -> all times)
+ `visible`  (optional) Figure visibility (default: "on")
+ `triadLineWidth`  (optional) Line width for triad/inertial flux lines (default: 2)
+ `forcingLineWidth`  (optional) Line width for forcing flux lines (default: 1.5)
+ `fluxTolerance`  (optional) Minimum absolute flux magnitude to display (default: 5e-2)
+ `keFractions`  (optional) Fractional ticks for KE/(KE+APE) axis (default: [.01,0.025,.1,.25,.5,.75,.9,0.975,.99])
+ `forcingFluxAttributes`  Array/struct to override forcing attributes (color, fancyName, alpha)

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency
  Produces one-dimensional summaries of energy fluxes (forcing and inertial)
  projected onto a radial (pseudo-)wavenumber axis, KE/(KE+APE) axis, and
  frequency axis. Computes temporal averages over the requested timeIndices,
  prepares per-flux metadata (colors, styles, alpha), and draws a compact
  2x2 multipanel figure that highlights triad and forcing contributions.
 
                      
