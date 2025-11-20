---
layout: default
title: logWavelengthAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 65
mathjax: true
---

#  logWavelengthAxis

Produce tick labels and positions for a log-scaled wavelength x-axis.


---

## Declaration
```matlab
 [labels, ticks] = logWavelengthAxis(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.num_ticks`  (optional) number of tick labels to generate (default: 6)
+ `options.roundToNearest`  (optional) round wavelength labels to this nearest integer (km) (default: 5)

## Returns
+ `labels`  cell array of label strings (km)
+ `ticks`  numeric vector of tick positions (units match the x-axis used in spectral plots)

## Discussion

  Returns formatted label strings (pseudo-wavelength in km) and numeric
  tick positions suitable for use with xticks/xticklabels on spectral
  figures. The number of ticks and rounding of the displayed wavelength
  values are controlled via options.
 
              
