---
layout: default
title: overlayGeostrophicKineticPotentialRatioContours
parent: WVDiagnostics
grand_parent: Classes
nav_order: 71
mathjax: true
---

#  overlayGeostrophicKineticPotentialRatioContours

Overlay contours of geostrophic kinetic-to-potential energy ratio.


---

## Declaration
```matlab
 overlayGeostrophicKineticPotentialRatioContours(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.ratios`  vector of ratio levels to contour (default: [.25 .5 .75])
+ `options.lineWidth`  contour line width (default: 1)
+ `options.textColor`  color for contour labels (default: [0.7 0.7 0.7])
+ `options.labelSpacing`  label spacing passed to clabel (default: 1000)

## Discussion

  Adds contours showing the local ratio KE_g/(KE_g+PE_g) on the
  current axes. Useful to indicate balanced vs unbalanced regions on
  energy-spectrum plots.
 
                - Returns: None. Draws contours on the current axes.
