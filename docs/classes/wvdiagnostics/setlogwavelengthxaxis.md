---
layout: default
title: setLogWavelengthXAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 116
mathjax: true
---

#  setLogWavelengthXAxis

Configure axis tick labels for a log-scaled wavelength x-axis.


---

## Declaration
```matlab
 setLogWavelengthXAxis(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.num_ticks`  (optional) Number of tick labels to generate (default: 6)
+ `options.roundToNearest`  (optional) Round wavelength labels to this nearest integer (km) (default: 5)

## Discussion

  Converts internal radial-wavenumber tick values into human-readable
  pseudo-wavelength labels (kilometers) and applies them to the current
  axes. Useful for spectral plots where the x-axis uses a log-scaled
  wavenumber coordinate but labels should show wavelengths in km.
 
            - Returns: None. Sets xticks and xticklabels on current axes.
