---
layout: default
title: plotEnergySpectrum
parent: WVDiagnostics
grand_parent: Classes
nav_order: 78
mathjax: true
---

#  plotEnergySpectrum

Plot the wave/geostrophic energy spectra at a given time.


---

## Declaration
```matlab
 [fig, spectralSlopes] = plotEnergySpectrum(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `iTime`  time index in model output file
+ `visible`  (optional) figure visibility (default: "on")

## Returns
+ `fig`  handle to the generated figure
+ `spectralSlopes`  struct containing spectral slope fits
  (`IOIGW_kR_slope`, `A0_kR_slope`, `IOIGW_j_slope`, `A0_j_slope`,
  `IOIGW_jWavenumber_slope`, `A0_jWavenumber_slope`). Fit lines are shown
  on the 1D spectra when this output is requested.

## Discussion

  Plot the wave/geostrophic energy spectra at a given time
  Makes a nice multiplanel plot of the wave and geostrophic spectra at a
  given time.
 
            
