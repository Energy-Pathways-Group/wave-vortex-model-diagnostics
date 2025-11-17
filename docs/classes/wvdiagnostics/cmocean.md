---
layout: default
title: cmocean
parent: WVDiagnostics
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  cmocean

returns perceptually-uniform colormaps created by Kristen Thyng.


---

## Declaration
```matlab
 cmap = cmocean(ColormapName,varargin)
```
## Parameters
+ `ColormapName`  input argument `ColormapName`
+ `varargin`  input argument `varargin`

## Returns
+ `cmap`  output value `cmap`

## Discussion

  cmocean returns perceptually-uniform colormaps created by Kristen Thyng.
  % Syntax
  cmocean
  cmap = cmocean('ColormapName')
  cmap = cmocean('-ColormapName')
  cmap = cmocean(...,NLevels)
  cmap = cmocean(...,'pivot',PivotValue)
  cmap = cmocean(...,'negative')
  cmocean(...)
  % Description
  cmocean without any inputs displays the options for colormaps.
  cmap = cmocean('ColormapName') returns a 256x3 colormap. ColormapName can be any of
  of the following:
  SEQUENTIAL:                DIVERGING:
  'thermal'                  'balance'
  'haline'                   'delta'
  'solar'                    'curl'
  'ice'                      'diff'
  'gray'                     'tarn'
  'oxy'
  'deep'                     CONSTANT LIGHTNESS:
  'dense'                    'phase'
  'algae'
  'matter'                   OTHER:
  'turbid'                   'topo'
  'speed'
  'amp'
  'tempo'
  'rain'
  cmap = cmocean('-ColormapName') a minus sign preceeding any ColormapName flips the
  order of the colormap.
  cmap = cmocean(...,NLevels) specifies a number of levels in the colormap.  Default
  value is 256.
  cmap = cmocean(...,'pivot',PivotValue) centers a diverging colormap such that white
  corresponds to a given value and maximum extents are set using current caxis limits.
  If no PivotValue is set, 0 is assumed. Early versions of this function used 'zero'
  as the syntax for 'pivot',0 and the old syntax is still supported.
  cmap = cmocean(...,'negative') inverts the lightness profile of the colormap. This can be
  useful particularly for divergent colormaps if the default white point of divergence
  gets lost in a white background.
  cmocean(...) without any outputs sets the current colormap to the current axes.
  % Examples
  Using this sample plot:
  imagesc(peaks(1000)+1)
  colorbar
  Set the colormap to 'algae':
  cmocean('algae')
  Same as above, but with an inverted algae colormap:
  cmocean('-algae')
  Set the colormap to a 12-level 'solar':
  cmocean('solar',12)
  Get the RGB values of a 5-level thermal colormap:
  RGB = cmocean('thermal',5)
  Some of those values are below zero and others are above. If this dataset represents
  anomalies, perhaps a diverging colormap is more appropriate:
  cmocean('balance')
  It's unlikely that 1.7776 is an interesting value about which the data values
  diverge.  If you want to center the colormap on zero using the current color
  axis limits, simply include the 'pivot' option:
  cmocean('balance','pivot',0)
  % Author Info
  This function was written by Chad A. Greene of the Institute for Geophysics at the
  University of Texas at Austin (UTIG), June 2016, using colormaps created by Kristen
  Thyng of Texas A&M University, Department of Oceanography. More information on the
  cmocean project can be found at http://matplotlib.org/cmocean/.
  % Citing this colormap:
  If you find an occasion to cite these colormaps for any reason, or if you just want
  some nice beach reading, check out the following paper from the journal Oceanography:
  Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True
  colors of oceanography: Guidelines for effective and accurate colormap selection.
  Oceanography 29(3):9?13, http://dx.doi.org/10.5670/oceanog.2016.66.
  See also colormap and caxis.
 
          
