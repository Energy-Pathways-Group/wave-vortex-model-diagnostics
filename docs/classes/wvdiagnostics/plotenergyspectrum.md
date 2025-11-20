---
layout: default
title: plotEnergySpectrum
parent: WVDiagnostics
grand_parent: Classes
nav_order: 76
mathjax: true
---

#  plotEnergySpectrum

Plot the wave/geostrophic energy spectra at a given time.


---

## Declaration
```matlab
 fig = plotEnergySpectrum(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `iTime`  time index in model output file
+ `visible`  (optional) figure visibility (default: "on")

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot the wave/geostrophic energy spectra at a given time
  Makes a nice multiplanel plot of the wave and geostrophic spectra at a
  given time.
 
            
