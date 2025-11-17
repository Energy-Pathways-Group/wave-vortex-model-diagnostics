---
layout: default
title: plotMooringRotarySpectrum
parent: WVDiagnostics
grand_parent: Classes
nav_order: 88
mathjax: true
---

#  plotMooringRotarySpectrum

Plot rotary spectra from mooring velocity time series.


---

## Declaration
```matlab
 fig = plotMooringRotarySpectrum(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `visible`  (optional) Figure visibility (default: "on")
+ `title`  (optional) Overall figure title (default: "total velocity"; set to "none" to suppress)
+ `shouldShowLegend`  (optional) (optional, logical) Show legend (default: true)
+ `shouldShowSpectralTitles`  (optional) (optional, logical) Show per-panel titles (default: true)

## Returns
+ `fig`  Handle to the generated figure

## Discussion

  Plot rotary spectra from mooring velocity time series.
  This function assuminges mooring horizontal velocity time series are saved to file.
  This function uses mspec from jLab.
  Compute and plot negative and positive rotary spectra from mooring
  horizontal velocity records stored in the WVModel output. Reads mooring
  time and velocity variables from self.wvfile, computes multitaper rotary
  spectra for each depth, averages over mooring IDs, and draws a two-panel
  log-log figure of negative and positive rotary spectra. Annotates plots
  with inertial frequency (f), M2 tidal frequency, f+M2, and min/max buoyancy
  frequency (N). Legend, panel titles, and overall title are optional.
 
                
