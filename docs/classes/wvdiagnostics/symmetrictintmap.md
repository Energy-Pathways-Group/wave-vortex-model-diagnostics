---
layout: default
title: symmetricTintMap
parent: WVDiagnostics
grand_parent: Classes
nav_order: 129
mathjax: true
---

#  symmetricTintMap

Colormap going from color -> tinted white -> color


---

## Discussion

  This function is used by plotPoissonFlowOverContours
 
    c   : 1x3 RGB color in [0 1]
    N   : total number of levels (even recommended)
    tintStrength : fraction of c mixed into white at the center
 
  Example:
    c = [0.2 0.6 0.8];
    colormap(symmetricTintMap(c,256,0.05)); colorbar
