---
layout: default
title: plotPoissonFlowOverContours
parent: WVDiagnostics
grand_parent: Classes
nav_order: 89
mathjax: true
---

#  plotPoissonFlowOverContours

Plot Poisson Flow Over Contours.


---

## Declaration
```matlab
 fig = plotPoissonFlowOverContours(wvd,options)
```
## Parameters
+ `wvd`  WVDiagnostics object
+ `visible`  (optional) input argument `visible` (default: "on")
+ `inertialFlux`  input argument `inertialFlux`
+ `vectorDensityLinearTransitionWavenumber`  (optional) input argument `vectorDensityLinearTransitionWavenumber` (default: Inf)
+ `forcingFlux`  input argument `forcingFlux`
+ `wavelengths`  (optional) input argument `wavelengths` (default: [1,2,5,10,20,50,100,200,500])
+ `wavelengthColor`  (optional) input argument `wavelengthColor` (default: [.5,.5,.5])
+ `addFrequencyContours`  (optional) input argument `addFrequencyContours` (default: false)
+ `frequencies`  (optional) input argument `frequencies` (default: [1.01 1.05 1.2 2 4 8 16])
+ `frequencyColor`  (optional) input argument `frequencyColor` (default: [.7,.7,.7])
+ `addKEPEContours`  (optional) input argument `addKEPEContours` (default: false)
+ `keFractions`  (optional) input argument `keFractions` (default: [.01,.1,.25,.5,.75,.9,.99])
+ `keFractionColor`  (optional) input argument `keFractionColor` (default: [.7,.7,.7])
+ `labelSpacing`  (optional) input argument `labelSpacing` (default: 1000)
+ `lineWidth`  (optional) input argument `lineWidth` (default: 1)
+ `kmax`  (optional) input argument `kmax` (default: Inf)
+ `jmax`  (optional) input argument `jmax` (default: Inf)
+ `quiverScale`  input argument `quiverScale`
+ `figureHandle`  input argument `figureHandle`
+ `nLevels`  (optional) input argument `nLevels` (default: 10)

## Returns
+ `fig`  Figure handle for the generated plot

## Discussion

  plotPoissonFlowOverContours is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
 
                                              
