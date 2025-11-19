---
layout: default
title: plotFluidStateMultipanel
parent: WVDiagnostics
grand_parent: Classes
nav_order: 83
mathjax: true
---

#  plotFluidStateMultipanel

Plot multipanel summary of fluid state and spectra


---

## Declaration
```matlab
 fig = plotFluidStateMultipanel(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `visible`  (optional) input argument `visible` (default: "on")
+ `iTime`  time-related parameter `iTime`
+ `title`  input argument `title`
+ `shouldShowEnergySpectra`  (optional) input argument `shouldShowEnergySpectra` (default: true)
+ `shouldShowTotalFields`  (optional) input argument `shouldShowTotalFields` (default: false)
+ `figureHandle`  input argument `figureHandle`
+ `wavelengths`  (optional) input argument `wavelengths` (default: [1,2,5,10,20,50,100,200,500])
+ `wavelengthColor`  (optional) input argument `wavelengthColor` (default: [.5,.5,.5])
+ `frequencies`  (optional) input argument `frequencies` (default: [1.01 1.05 1.2 2 4 8 16])
+ `frequencyColor`  (optional) input argument `frequencyColor` (default: [.7,.7,.7])
+ `keFractions`  (optional) input argument `keFractions` (default: [.01,.1,.25,.5,.75,.9,.99])
+ `keFractionColor`  (optional) input argument `keFractionColor` (default: [.7,.7,.7])
+ `labelSpacing`  (optional) input argument `labelSpacing` (default: 1000)
+ `lineWidth`  (optional) input argument `lineWidth` (default: 1)

## Returns
+ `fig`  Figure handle for the generated plot

## Discussion

  Create a compact multipanel figure showing horizontal (x-y) maps and
  vertical (x-z) sections of vertical vorticity for the total flow, the
  wave component, and the geostrophic component at a specified model time.
  Optionally includes log-energy spectra with KE/PE and frequency/wavelength
  contours. Axes are annotated in kilometers and depth in kilometers; color
  limits and annotation ticks are handled internally.
 
                                    
