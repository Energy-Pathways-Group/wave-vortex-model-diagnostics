---
layout: default
title: overlayGeostrophicKineticPotentialFractionContours
parent: WVDiagnostics
grand_parent: Classes
nav_order: 70
mathjax: true
---

#  overlayGeostrophicKineticPotentialFractionContours

Overlay contour where geostrophic kinetic fraction equals a given value.


---

## Declaration
```matlab
 overlayGeostrophicKineticPotentialFractionContours(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.fraction`  fractional contour level to draw (default: 0.5)
+ `options.textColor`  color for contour labels (default: [0.7 0.7 0.7])
+ `options.labelSpacing`  label spacing passed to clabel (default: 1000)

## Discussion

  Draws the contour line where KE_g / (KE_g + PE_g) equals a chosen
  fractional value (default 0.5). Useful to mark equipartition lines
  on spectral diagrams.
 
              - Returns: None. Draws contour on the current axes.
