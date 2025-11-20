---
layout: default
title: showRossbyRadiusYAxis
parent: WVDiagnostics
grand_parent: Classes
nav_order: 119
mathjax: true
---

#  showRossbyRadiusYAxis

Annotate the axes with a Rossby radius y-axis label.


---

## Declaration
```matlab
 showRossbyRadiusYAxis(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `options.textColor`  color for the annotation text (default: 'k')
+ `options.xFrac`  fractional x position in axes coordinates (default: 0.02)
+ `options.yFrac`  fractional y position in axes coordinates (default: 0.5)

## Discussion

  Places a textual annotation (L_r) on the current axes indicating the
  Rossby radius scale in kilometers. Positioning and color can be
  adjusted via options.
 
              - Returns: handle to the created text object
