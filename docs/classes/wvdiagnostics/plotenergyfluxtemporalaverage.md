---
layout: default
title: plotEnergyFluxTemporalAverage
parent: WVDiagnostics
grand_parent: Classes
nav_order: 74
mathjax: true
---

#  plotEnergyFluxTemporalAverage

Plot temporally averaged energy flux diagnostics.


---

## Declaration
```matlab
 fig = plotEnergyFluxTemporalAverage(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `approximation`  (optional) {'quadratic','exact'} approximation to use (default: 'quadratic')
+ `energyReservoir`  (optional) EnergyReservoir to plot (default: EnergyReservoir.total)
+ `triadComponents`  (optional) TriadFlowComponent vector for inertial flux selection (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
+ `showForcingFluxes`  (optional) (optional, logical) Include forcing fluxes when using quadratic approximation (default: true)
+ `timeIndices`  (optional) Time indices to average over (default: Inf -> all times)
+ `axes`  (optional) Plot projection, one of {'jk','j','jWavenumber','k','k-pseudo-isotropic','omega'} (default: 'jk')
+ `filter`  (optional) Function handle applied to plotted data (default: @(v) v)
+ `shouldOverlayWaveFrequencies`  (optional) (optional, logical) Overlay wave-frequency contours on 'jk' axes (default: false)
+ `shouldOverlayGeostrophicKineticPotentialRatioContours`  (optional) (optional, logical) Overlay geostrophic KE/PE fraction contours (default: true)
+ `colormap`  (optional) Colormap to use for 'jk' axes (default: WVDiagnostics.crameri('-bam'))
+ `visible`  (optional) Figure visibility (default: "on")
+ `overSaturationFactor`  (optional) Scalar or two-element colormap limits (default: 10)
+ `fluxGroups`  (optional) Cell array of index groups to combine for plotting (default: [])
+ `simpleName`  (optional) Cell array of simple names for fluxGroups (default: [])

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot temporally averaged energy flux diagnostics.
  Compute and plot temporally averaged energy flux diagnostics (exact or
  quadratic approximation) on a chosen axes projection. Supports grouping
  fluxes, overlaying frequency/ratio contours, and customizing colormap,
  saturation, and figure visibility.
 
                                    
