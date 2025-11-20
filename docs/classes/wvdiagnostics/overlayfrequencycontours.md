---
layout: default
title: overlayFrequencyContours
parent: WVDiagnostics
grand_parent: Classes
nav_order: 68
mathjax: true
---

#  overlayFrequencyContours

Overlay contours of nondimensional frequency on the current axes.


---

## Declaration
```matlab
 overlayFrequencyContours(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.frequencies`  vector of nondimensional frequency values (omega/f) to draw (default: [])
+ `options.lineWidth`  contour line width (default: 1)
+ `options.textColor`  color for contour labels (default: [0.7 0.7 0.7])
+ `options.labelSpacing`  label spacing passed to clabel (default: 1000)

## Discussion

  Draws contours of the model dispersion frequency (omega/f) converted
  to the plotting coordinates and optionally labels them. Intended for
  spectrum figures to show dispersion curves.
 
                - Returns: None. Draws contours on the current axes.
