---
layout: default
title: WVDiagnostics
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 1
---

#  WVDiagnostics

Produces diagnostics and figures from WVModel output


---

## Overview
    This is a collection of diagnostic tools for analyzing model output
 
                                                      


## Topics
+ Diagnostics Generation
  + [`create1DMirrorFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/create1dmirrorfluxes.html) Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
  + [`create2DMirrorFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/create2dmirrorfluxes.html) Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
  + [`createDiagnosticsFile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/creatediagnosticsfile.html) Create a new diagnostics file and compute diagnostics from WVModel output.
  + [`createReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/createreservoirgroup.html) Create a reservoir group in the diagnostics NetCDF and populate it.
+ Figures
  + Model Snapshot
    + [`plotEnergySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyspectrum.html) Plot the wave/geostrophic energy spectra at a given time.
    + [`plotEnstrophySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyspectrum.html) Plot the enstrophy spectrum at a given time.
    + [`plotFluidDecompositionMultipanel`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotfluiddecompositionmultipanel.html) Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
    + [`plotFluidStateMultipanel`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotfluidstatemultipanel.html) Plot multipanel summary of fluid state and spectra
    + [`plotPotentialEnergySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpotentialenergyspectrum.html) Plot the potential energy spectrum
  + Energy
    + [`plotEnergyFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxovertime.html) Plot energy fluxes as a function of time.
    + [`plotEnergyFluxTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxtemporalaverage.html) Plot temporally averaged energy flux diagnostics.
    + [`plotEnergyFluxes1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxes1d.html) Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency.
    + [`plotEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyovertime.html) Plot energy for each reservoir over time.
    + [`plotEnergyTriadFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergytriadfluxovertime.html) Plot inertial flux for each reservoir over time
    + [`plotSourcesSinksForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotsourcessinksforreservoirgroup.html) Plot sources, sinks, and reservoirs diagram for a reservoir group.
    + [`plotSourcesSinksReservoirsDiagram`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotsourcessinksreservoirsdiagram.html) Plot sources, sinks, and reservoirs diagram.
    + [`showDampingFluxVsPseudolength`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/showdampingfluxvspseudolength.html) Show Damping Flux Vs Pseudolength.
  + Potential Enstrophy
    + [`plotEnstrophyFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyfluxovertime.html) Plot enstrophy fluxes for reservoirs as a function of time.
    + [`plotEnstrophyFluxTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyfluxtemporalaverage.html) Plot temporally averaged enstrophy flux diagnostics.
    + [`plotEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyovertime.html) Plot enstrophy over time
    + [`plotEnstrophyTriadFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophytriadfluxovertime.html) Plot enstrophy triad fluxes over time.
  + Ancillary
    + [`plotMooringRotarySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotmooringrotaryspectrum.html) Plot rotary spectra from mooring velocity time series.
  + Auxiliary functions
    + [`logWavelengthAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/logwavelengthaxis.html) Produce tick labels and positions for a log-scaled wavelength x-axis.
    + [`overlayFrequencyContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlayfrequencycontours.html) Overlay contours of nondimensional frequency on the current axes.
    + [`overlayGeostrophicKineticPotentialFractionContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlaygeostrophickineticpotentialfractioncontours.html) Overlay contour where geostrophic kinetic fraction equals a given value.
    + [`overlayGeostrophicKineticPotentialRatioContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlaygeostrophickineticpotentialratiocontours.html) Overlay contours of geostrophic kinetic-to-potential energy ratio.
    + [`setLogWavelengthXAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/setlogwavelengthxaxis.html) Configure axis tick labels for a log-scaled wavelength x-axis.
    + [`showRossbyRadiusYAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/showrossbyradiusyaxis.html) Annotate the axes with a Rossby radius y-axis label.
  + Diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`plotPoissonFlowOverContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpoissonflowovercontours.html) Plot Poisson Flow Over Contours.
+ Summaries
  + [`createEnstrophyFluxSummaryTable`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/createenstrophyfluxsummarytable.html) Create a LaTeX table summarizing enstrophy source/sink fluxes.
  + [`summarizeSourcesSinksReservoirs`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/summarizesourcessinksreservoirs.html) Summarize Sources Sinks Reservoirs.
+ Diagnostics
  + Energy
    + Time series, [t 1]
      + [`exactEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyovertime.html) Exact Energy Over Time.
      + [`quadraticEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyovertime.html) Compute energy for each reservoir over time.
  + Energy Fluxes
    + General, [j kRadial t]
      + [`exactEnergyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxes.html) Return the exact energy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnergyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxes.html) Return the energy flux from the forcing terms.
      + [`quadraticEnergyTriadFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxes.html) Return the energy flux from the inertial terms, specified as triad components.
    + Temporal averages, [j kRadial]
      + [`exactEnergyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxestemporalaverage.html) Temporally averaged exact energy fluxes.
      + [`quadraticEnergyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxestemporalaverage.html) Compute temporally averaged forcing fluxes.
      + [`quadraticEnergyTriadFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxestemporalaverage.html) Computes the temporally averaged quadratic triad fluxes.
    + Time series, [t 1]
      + [`exactEnergyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxesovertime.html) Exact energy fluxes over time.
      + [`quadraticEnergyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxesovertime.html) Compute forcing fluxes over time.
    + Spatial-temporal averages, [1 1]
      + [`exactEnergyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact forcing fluxes.
      + [`quadraticEnergyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of forcing fluxes.
      + [`quadraticEnergyTriadFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxesspatialtemporalaverage.html) Compute spatial-temporal average of inertial fluxes.
    + Temporal averages, 1D axes
      + [`quadraticEnergyMirrorTriadFluxes1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d.html) Quadratic Energy Mirror Triad Fluxes1 D.
      + [`quadraticEnergyMirrorTriadFluxes1D_kepe`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_kepe.html) Quadratic Energy Mirror Triad Fluxes1 D kepe.
      + [`quadraticEnergyMirrorTriadFluxes1D_omega`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_omega.html) Quadratic Energy Mirror Triad Fluxes1 D omega.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_kepe.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D kepe.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_omega.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D omega.
    + Temporal averages, 2D axes [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyMirrorTriadFluxes2D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes2d.html) Quadratic Energy Mirror Triad Fluxes2 D.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage2D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage2d.html) outputGrid determines whether or not the fluxes get downsampled to the.
  + Potential Enstrophy
    + Time series, [t 1]
      + [`exactEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyovertime.html) Exact Enstrophy Over Time.
      + [`quadraticEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyovertime.html) Quadratic Enstrophy Over Time.
  + Potential Enstrophy Fluxes
    + General, [j kRadial t]
      + [`exactEnstrophyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxes.html) Return the available potential enstrophy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnstrophyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxes.html) Return the enstrophy flux from the forcing terms.
      + [`quadraticEnstrophyTriadFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophytriadfluxes.html) Return the enstrophy flux from the forcing terms.
    + Temporal averages, [j kRadial]
      + [`exactEnstrophyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
      + [`quadraticEnstrophyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
    + Time series, [t 1]
      + [`exactEnstrophyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxesovertime.html) Compute exact enstrophy fluxes over time.
      + [`quadraticEnergyTriadFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxesovertime.html) Compute inertial fluxes over time.
      + [`quadraticEnstrophyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxesovertime.html) Compute enstrophy fluxes over time.
    + Spatial-temporal averages, [1 1]
      + [`exactEnstrophyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact enstrophy fluxes.
      + [`quadraticEnstrophyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the qgpv enstrophy fluxes.
  + Flux diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`RescalePoissonFlowFluxForLogSpace`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/rescalepoissonflowfluxforlogspace.html) Rescale Poisson Flow Flux For Log Space.
  + General — Misc
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`cmocean`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cmocean.html) returns perceptually-uniform colormaps created by Kristen Thyng.
      + [`geostrophicGeostrophicWaveEnergy`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/geostrophicgeostrophicwaveenergy.html) Note that.
      + [`waveWaveGeostrophicEnergy`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wavewavegeostrophicenergy.html) Note that.
      + [`waveWaveGeostrophicEnergyForMode`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wavewavegeostrophicenergyformode.html) Note that.
  + Enstrophy fluxes — Triad interactions
    + Fluxes over time, [t 1]
      + [`quadraticEnstrophyTriadFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophytriadfluxesovertime.html) Compute enstrophy inertial (aka, triad) fluxes over time.
+ Transforms
  + Spectral
    + General
+ Utilities
  + Colormaps
    + Crameri
  + Sparse matrices
    + Axis binning
  + Colormaps — Crameri
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`crameri`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/crameri.html) returns perceptually-uniform scientific colormaps created.
  + Sparse matrices — Axis binning
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`sparseJKAxisBinMatrices`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparsejkaxisbinmatrices.html) Sparse JKAxis Bin Matrices.
      + [`sparseJWavenumberAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparsejwavenumberaxis.html) Sparse JWavenumber Axis.
      + [`sparseKRadialAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparsekradialaxis.html) Sparse KRadial Axis.
      + [`sparseKePeAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparsekepeaxis.html) Sparse Ke Pe Axis.
      + [`sparseOmegaAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparseomegaaxis.html) Sparse Omega Axis.
      + [`sparsePseudoRadialAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sparsepseudoradialaxis.html) Sparse Pseudo Radial Axis.
+ Internal
  + Support functions for createReservoirGroup
    + [`addTriadFluxesForReservoirGroupAtTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/addtriadfluxesforreservoirgroupattime.html) Add Triad Fluxes For Reservoir Group At Time.
    + [`filterEnergyForSourcesSinksReservoirs`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/filterenergyforsourcessinksreservoirs.html) This function returns values assuming three reservoirs: geo, wave, and.
    + [`fluxesForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/fluxesforreservoirgroup.html) Fluxes For Reservoir Group.
    + [`variablesForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/variablesforreservoirgroup.html) Variables For Reservoir Group.
+ Transformations
  + Cosine and Sine
    + [`CosineTransformBack`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cosinetransformback.html) Fast Discrete Inverse Cosine Transform.
    + [`CosineTransformForward`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cosinetransformforward.html) Fast Discrete Cosine Transform (DCT-I).
    + [`DCT1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dct1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
    + [`DCT2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dct2.html) CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix.
    + [`DST1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dst1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
    + [`DST2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dst2.html) SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix.
    + [`SineTransformBack`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sinetransformback.html) Fast Discrete Inverse Sine Transform.
    + [`iDCT1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idct1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
    + [`iDCT2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idct2.html) InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix.
    + [`iDST1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idst1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
    + [`iDST2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idst2.html) InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix.
    + General
      + [`SineTransformForward`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sinetransformforward.html) Fast Discrete Sine Transform (DST-I).
  + Axes
    + [`transformToKePeAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtokepeaxis.html) Transform To Ke Pe Axis.
    + [`transformToOmegaAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtoomegaaxis.html) transforms in the from (j,kRadial) to omegaAxis.
    + [`transformToPseudoRadialWavenumber`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumber.html) transforms in the from (j,kRadial) to kPseudoRadial.
    + [`transformToPseudoRadialWavenumberA0`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumbera0.html) transforms in the from (j,kRadial) to kPseudoRadial.
    + [`transformToPseudoRadialWavenumberApm`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumberapm.html) transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial.
+ Constructor
  + [`WVDiagnostics`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvdiagnostics.html) Initializes the WVDiagnostics object, loads the wave-vortex transform and diagnostics files.
+ Computing Spectra
  + [`crossSpectrumWithFgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/crossspectrumwithfgtransform.html) Cross Spectrum With Fg Transform.
  + [`crossSpectrumWithGgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/crossspectrumwithggtransform.html) Cross Spectrum With Gg Transform.
  + [`spectrumWithFgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/spectrumwithfgtransform.html) Spectrum With Fg Transform.
  + [`spectrumWithGgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/spectrumwithggtransform.html) Spectrum With Gg Transform.
+ Units
  + [`escale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/escale.html) energy scaling (numeric value)
  + [`escale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/escale_units.html) energy units (string value)
  + [`flux_scale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/flux_scale.html) energy flux scaling (numeric value)
  + [`flux_scale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/flux_scale_units.html) energy flux units (string value)
  + [`setEnergyUnits`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/setenergyunits.html) Set the time and energy scaling and units for plotting and output.
  + [`tscale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/tscale.html) time axis scaling (numeric value)
  + [`tscale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/tscale_units.html) time axis units (string value)
  + [`z_flux_scale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/z_flux_scale.html) potential enstrophy flux scaling (numeric value)
  + [`z_flux_scale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/z_flux_scale_units.html) potential enstrophy flux units (string value)
  + [`zscale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/zscale.html) potential enstrophy scaling (numeric value)
  + [`zscale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/zscale_units.html) potential enstrophy units (string value)
+ Dependent property getter
  + [`t_diag`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/t_diag.html) Get time vector from the diagnostics file
+ Other
  + [`Lr2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/lr2.html) 
  + [`Lr2_pm`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/lr2_pm.html) 
  + [`PoissonFlowFromFlux`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromflux.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxDCTI`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxdcti.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxType1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxtype1.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxWithAxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxwithaxes.html) We will treat the first dimension as `x' and the second
  + [`areEnergyReservoirsComplete`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/areenergyreservoirscomplete.html) 
  + [`areTriadComponentsComplete`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/aretriadcomponentscomplete.html) 
  + [`diagfile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagfile.html) NetCDFFile instance for the diagnostics file
  + [`diagnosticsHasExplicitAntialiasing`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagnosticshasexplicitantialiasing.html) 
  + [`diagpath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagpath.html) path to the diagnostics file
  + [`fancyNameForName`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/fancynameforname.html) 
  + [`filterFluxesForReservoir`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/filterfluxesforreservoir.html) 
  + [`forcingNames`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/forcingnames.html) 
  + [`geo_hke_jk`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/geo_hke_jk.html) 
  + [`geo_pe_jk`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/geo_pe_jk.html) 
  + [`iTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/itime.html) 
  + [`iTimeChanged`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/itimechanged.html) 
  + [`j`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/j.html) 
  + [`jWavenumber`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/jwavenumber.html) 
  + [`kPseudoRadial`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/kpseudoradial.html) 
  + [`kRadial`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/kradial.html) 
  + [`kePeAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/kepeaxis.html) 
  + [`omegaAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/omegaaxis.html) 
  + [`omega_jk`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/omega_jk.html) 
  + [`plotPoissonFlowOverPcolor`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpoissonflowoverpcolor.html) 
  + [`pseudoRadialBinning`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/pseudoradialbinning.html) choose an algorithm for pseudoRadialBinning
  + [`symmetricTintMap`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/symmetrictintmap.html) Colormap going from color -> tinted white -> color
  + [`t_wv`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/t_wv.html) 
  + [`version`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/version.html) Locate folder containing this class file
  + [`wvaapath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvaapath.html) path to a WVTransform with explicit anti-aliasing
  + [`wvfile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvfile.html) NetCDFFile instance for the WaveVortexModel output
  + [`wvpath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvpath.html) path to the WaveVortexModel output
  + [`wvt`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvt.html) WVTransform instance, set to the current iTime
  + [`wvt_aa`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvt_aa.html) 


---