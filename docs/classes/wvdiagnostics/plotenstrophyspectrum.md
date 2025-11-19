---
layout: default
title: plotEnstrophySpectrum
parent: WVDiagnostics
grand_parent: Classes
nav_order: 80
mathjax: true
---

#  plotEnstrophySpectrum

Plot the enstrophy spectrum at a given time.


---

## Declaration
```matlab
 fig = plotEnstrophySpectrum(self,options)
```
## Parameters
+ `self`  WVDiagnostics object
+ `iTime`  time index in model output file (default: self.iTime)
+ `visible`  (optional) Figure visibility (default: "on")
+ `title`  Figure title (default: 'Enstrophy Spectrum')
+ `clim`  (optional) Color limits for log10 pcolor panels (default: [-13 -7])
+ `figureHandle`  Existing figure handle to draw into (default: create new figure)
+ `style`  (optional) Plot style, one of "four-panel" or "kPseudoRadial" (default: "four-panel")

## Returns
+ `fig`  handle to the generated figure

## Discussion

  Plot the enstrophy spectrum at a given time
  Produces a multiplanel figure of geostrophic and apv enstrophy diagnostics at a
  specified time index. Shows available potential enstrophy (APV), the
  quadratic (QGPV) approximation, their difference, plus modal and radial
  spectra. Supports a 2x2 "four-panel" layout or a single-panel
  pseudo-radial projection ("kPseudoRadial").
 
                    
