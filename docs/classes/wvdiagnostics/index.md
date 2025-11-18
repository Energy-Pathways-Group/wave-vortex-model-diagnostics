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
 
  Topics overview:
    - Configuration — Reservoirs — Grouping
    - Diagnostics — Energy fluxes — General — Fluxes over time, [t 1]
    - Diagnostics — Energy fluxes — General — Fluxes, [j kRadial t]
    - Diagnostics — Energy fluxes — Time/space averages — Flux averages, scalar [1 1]
    - Diagnostics — Energy fluxes — Time/space averages — Fluxes in space, [j kRadial]
    - Diagnostics — Energy fluxes — Triad interactions — Flux averages, scalar [1 1]
    - Diagnostics — Energy fluxes — Triad interactions — Fluxes in space, [j kRadial]
    - Diagnostics — Energy fluxes — Triad interactions — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — Energy fluxes — Triad interactions — Fluxes over time, [t 1]
    - Diagnostics — Energy fluxes — Triad interactions — Fluxes, [j kRadial t]
    - Diagnostics — Energy fluxes — Triad interactions — Mirror pairs — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — Energy fluxes — Triad interactions — Mirror pairs — Fluxes in space, [sparseKRadialAxis 1]
    - Diagnostics — Energy fluxes — Triad interactions — Primary — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — Energy fluxes — Triad interactions — Primary — Fluxes in space, [sparseKRadialAxis 1]
    - Diagnostics — Energy — Time series
    - Diagnostics — Enstrophy fluxes — General
    - Diagnostics — Enstrophy fluxes — General — Fluxes over time, [t 1]
    - Diagnostics — Enstrophy fluxes — General — Fluxes, [j kRadial t]
    - Diagnostics — Enstrophy fluxes — Time/space averages — Flux averages, scalar [1 1]
    - Diagnostics — Enstrophy fluxes — Time/space averages — Fluxes in space, [j kRadial]
    - Diagnostics — Enstrophy fluxes — Triad interactions — Fluxes over time, [t 1]
    - Diagnostics — Enstrophy fluxes — Triad interactions — Fluxes, [j kRadial t]
    - Diagnostics — Enstrophy — Time series
    - Diagnostics — Flux diagnostics — General
    - Diagnostics — Flux diagnostics — General — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — Flux diagnostics — Triad interactions
    - Diagnostics — General — Misc
    - Diagnostics — General — Misc — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — General — Misc — Fluxes over time, [t 1]
    - Diagnostics — Spectra — Cross-spectra — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Diagnostics — Spectra — Potential energy — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Figures — Diagnostics — General
    - Figures — Diagnostics — General — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Figures — Spectra — Potential energy
    - Figures — Time series — Diagnostics
    - Transforms — Spectral — General
    - Transforms — Spectral — General — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Utilities — Colormaps — Crameri — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
    - Utilities — Sparse matrices — Axis binning — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]


## Topics
+ Transforms
  + Spectral — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`CosineTransformBack`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cosinetransformback.html) Fast Discrete Inverse Cosine Transform.
      + [`CosineTransformForward`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cosinetransformforward.html) Fast Discrete Cosine Transform (DCT-I).
      + [`SineTransformBack`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sinetransformback.html) Fast Discrete Inverse Sine Transform.
      + [`transformToKePeAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtokepeaxis.html) Transform To Ke Pe Axis.
      + [`transformToOmegaAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtoomegaaxis.html) transforms in the from (j,kRadial) to omegaAxis.
      + [`transformToPseudoRadialWavenumber`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumber.html) transforms in the from (j,kRadial) to kPseudoRadial.
      + [`transformToPseudoRadialWavenumberA0`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumbera0.html) transforms in the from (j,kRadial) to kPseudoRadial.
      + [`transformToPseudoRadialWavenumberApm`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/transformtopseudoradialwavenumberapm.html) transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial.
  + Spectral
    + General
      + [`SineTransformForward`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/sinetransformforward.html) Fast Discrete Sine Transform (DST-I).
+ Diagnostics
  + General — Misc
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`DCT1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dct1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`DCT2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dct2.html) CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix.
      + [`DST1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dst1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`DST2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/dst2.html) SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix.
      + [`cmocean`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/cmocean.html) returns perceptually-uniform colormaps created by Kristen Thyng.
      + [`geostrophicGeostrophicWaveEnergy`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/geostrophicgeostrophicwaveenergy.html) Note that.
      + [`iDCT1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idct1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`iDCT2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idct2.html) InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix.
      + [`iDST1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idst1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`iDST2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/idst2.html) InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix.
      + [`waveWaveGeostrophicEnergy`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wavewavegeostrophicenergy.html) Note that.
      + [`waveWaveGeostrophicEnergyForMode`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wavewavegeostrophicenergyformode.html) Note that.
    + Fluxes over time, [t 1]
      + [`quadraticEnergyMirrorTriadsUndamped`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadsundamped.html) Quadratic Energy Mirror Triads Undamped.
  + Flux diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`RescalePoissonFlowFluxForLogSpace`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/rescalepoissonflowfluxforlogspace.html) Rescale Poisson Flow Flux For Log Space.
  + Flux diagnostics
    + Triad interactions
      + [`addTriadFluxesForReservoirGroupAtTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/addtriadfluxesforreservoirgroupattime.html) Add Triad Fluxes For Reservoir Group At Time.
    + General
      + [`create1DMirrorFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/create1dmirrorfluxes.html) Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
      + [`create2DMirrorFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/create2dmirrorfluxes.html) Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
      + [`fluxesForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/fluxesforreservoirgroup.html) Fluxes For Reservoir Group.
      + [`showDampingFluxVsPseudolength`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/showdampingfluxvspseudolength.html) Show Damping Flux Vs Pseudolength.
  + General
    + Misc
      + [`createDiagnosticsFile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/creatediagnosticsfile.html) Create a new diagnostics file and compute diagnostics from WVModel output.
  + Enstrophy fluxes
    + General
      + [`createEnstrophyFluxSummaryTable`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/createenstrophyfluxsummarytable.html) Create a LaTeX table summarizing enstrophy source/sink fluxes.
  + Spectra — Cross-spectra
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`crossSpectrumWithFgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/crossspectrumwithfgtransform.html) Cross Spectrum With Fg Transform.
      + [`crossSpectrumWithGgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/crossspectrumwithggtransform.html) Cross Spectrum With Gg Transform.
  + Energy fluxes — General
    + Fluxes, [j kRadial t]
      + [`exactEnergyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxes.html) Return the exact energy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnergyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxes.html) Return the energy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`exactEnergyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxesovertime.html) Exact energy fluxes over time.
      + [`quadraticEnergyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxesovertime.html) Compute forcing fluxes over time.
  + Energy fluxes — Time/space averages
    + Flux averages, scalar [1 1]
      + [`exactEnergyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact forcing fluxes.
      + [`quadraticEnergyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of forcing fluxes.
    + Fluxes in space, [j kRadial]
      + [`exactEnergyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyfluxestemporalaverage.html) Temporally averaged exact energy fluxes.
      + [`quadraticEnergyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyfluxestemporalaverage.html) Compute temporally averaged forcing fluxes.
  + Energy
    + Time series
      + [`exactEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenergyovertime.html) Exact Energy Over Time.
      + [`quadraticEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyovertime.html) Compute energy for each reservoir over time.
  + Enstrophy fluxes — General
    + Fluxes, [j kRadial t]
      + [`exactEnstrophyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxes.html) Return the available potential enstrophy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnstrophyFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxes.html) Return the enstrophy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`exactEnstrophyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxesovertime.html) Compute exact enstrophy fluxes over time.
      + [`quadraticEnstrophyFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxesovertime.html) Compute enstrophy fluxes over time.
  + Enstrophy fluxes — Time/space averages
    + Flux averages, scalar [1 1]
      + [`exactEnstrophyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact enstrophy fluxes.
      + [`quadraticEnstrophyFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the qgpv enstrophy fluxes.
    + Fluxes in space, [j kRadial]
      + [`exactEnstrophyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
      + [`quadraticEnstrophyFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
  + Enstrophy
    + Time series
      + [`exactEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/exactenstrophyovertime.html) Exact Enstrophy Over Time.
      + [`quadraticEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophyovertime.html) Quadratic Enstrophy Over Time.
  + Energy fluxes — Triad interactions — Mirror pairs
    + Fluxes in space, [sparseKRadialAxis 1]
      + [`quadraticEnergyMirrorTriadFluxes1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d.html) Quadratic Energy Mirror Triad Fluxes1 D.
      + [`quadraticEnergyMirrorTriadFluxes1D_kepe`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_kepe.html) Quadratic Energy Mirror Triad Fluxes1 D kepe.
      + [`quadraticEnergyMirrorTriadFluxes1D_omega`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_omega.html) Quadratic Energy Mirror Triad Fluxes1 D omega.
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyMirrorTriadFluxes2D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergymirrortriadfluxes2d.html) Quadratic Energy Mirror Triad Fluxes2 D.
  + Energy fluxes — Triad interactions — Primary
    + Fluxes in space, [sparseKRadialAxis 1]
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_kepe.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D kepe.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_omega.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D omega.
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage2D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage2d.html) outputGrid determines whether or not the fluxes get downsampled to the.
  + Energy fluxes — Triad interactions
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyTriadFluxWWGWave`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxwwgwave.html) Quadratic Energy Triad Flux WWGWave.
    + Fluxes, [j kRadial t]
      + [`quadraticEnergyTriadFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxes.html) Return the energy flux from the inertial terms, specified as triad components.
    + Fluxes over time, [t 1]
      + [`quadraticEnergyTriadFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxesovertime.html) Compute inertial fluxes over time.
    + Flux averages, scalar [1 1]
      + [`quadraticEnergyTriadFluxesSpatialTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxesspatialtemporalaverage.html) Compute spatial-temporal average of inertial fluxes.
    + Fluxes in space, [j kRadial]
      + [`quadraticEnergyTriadFluxesTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenergytriadfluxestemporalaverage.html) Computes the temporally averaged quadratic triad fluxes.
  + Enstrophy fluxes — Triad interactions
    + Fluxes, [j kRadial t]
      + [`quadraticEnstrophyTriadFluxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophytriadfluxes.html) Return the enstrophy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`quadraticEnstrophyTriadFluxesOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/quadraticenstrophytriadfluxesovertime.html) Compute enstrophy inertial (aka, triad) fluxes over time.
  + Spectra — Potential energy
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`spectrumWithFgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/spectrumwithfgtransform.html) Spectrum With Fg Transform.
      + [`spectrumWithGgTransform`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/spectrumwithggtransform.html) Spectrum With Gg Transform.
+ Constructor
  + [`WVDiagnostics`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvdiagnostics.html) Initializes the WVDiagnostics object, loads the wave-vortex transform and diagnostics files.
+ Utilities
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
+ Configuration
  + [`setEnergyUnits`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/setenergyunits.html) Set the time and energy scaling and units for plotting and output.
  + Reservoirs
    + Grouping
      + [`createReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/createreservoirgroup.html) Create a reservoir group in the diagnostics NetCDF and populate it.
      + [`filterEnergyForSourcesSinksReservoirs`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/filterenergyforsourcessinksreservoirs.html) This function returns values assuming three reservoirs: geo, wave, and.
      + [`filterEnergyForSourcesSinksReservoirsOld`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/filterenergyforsourcessinksreservoirsold.html) This function returns values assuming three reservoirs: geo, wave, and.
      + [`summarizeSourcesSinksReservoirs`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/summarizesourcessinksreservoirs.html) Summarize Sources Sinks Reservoirs.
      + [`variablesForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/variablesforreservoirgroup.html) Variables For Reservoir Group.
+ Figures
  + Time series
    + Diagnostics
      + [`plotEnergyFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxovertime.html) Plot energy fluxes for reservoirs as a function of time.
      + [`plotEnergyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyovertime.html) Plot energy for each reservoir over time.
      + [`plotEnergyTriadFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergytriadfluxovertime.html) Plot inertial flux for each reservoir over time.
      + [`plotEnstrophyFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyfluxovertime.html) Plot enstrophy fluxes for reservoirs as a function of time.
      + [`plotEnstrophyOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyovertime.html) Plot enstrophy over time.
      + [`plotEnstrophyTriadFluxOverTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophytriadfluxovertime.html) Plot enstrophy triad fluxes over time.
  + Diagnostics
    + General
      + [`plotEnergyFluxTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxtemporalaverage.html) Plot temporally averaged energy flux diagnostics.
      + [`plotEnergyFluxes1D`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyfluxes1d.html) Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency.
      + [`plotEnstrophyFluxTemporalAverage`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyfluxtemporalaverage.html) Plot temporally averaged enstrophy flux diagnostics.
      + [`plotFluidDecompositionMultipanel`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotfluiddecompositionmultipanel.html) Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
      + [`plotFluidStateMultipanel`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotfluidstatemultipanel.html) Plot multipanel summary of fluid state and spectra
      + [`plotSourcesSinksForReservoirGroup`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotsourcessinksforreservoirgroup.html) Plot sources, sinks, and reservoirs diagram for a reservoir group.
      + [`plotSourcesSinksReservoirsDiagram`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotsourcessinksreservoirsdiagram.html) Plot sources, sinks, and reservoirs diagram.
      + [`plotSourcesSinksReservoirsDiagramWithClosureRegion`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotsourcessinksreservoirsdiagramwithclosureregion.html) Plot sources, sinks, and reservoirs diagram.
  + Spectra
    + Potential energy
      + [`plotEnergySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyspectrum.html) Plot the wave/geostrophic energy spectra at a given time.
      + [`plotEnergySpectrumOLD`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenergyspectrumold.html) Plot the wave/geostrophic energy spectra at a given time.
      + [`plotEnstrophySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyspectrum.html) Plot the enstrophy spectrum at a given time.
      + [`plotEnstrophySpectrumOLD`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotenstrophyspectrumold.html) Plot the geostrophic enstrophy spectrum at a given time.
      + [`plotMooringRotarySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotmooringrotaryspectrum.html) Plot rotary spectra from mooring velocity time series.
      + [`plotPotentialEnergySpectrum`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpotentialenergyspectrum.html) Plot energy for each reservoir over time.
  + Diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`plotPoissonFlowOverContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpoissonflowovercontours.html) Plot Poisson Flow Over Contours.
+ Other
  + [`Lr2`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/lr2.html) 
  + [`Lr2_pm`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/lr2_pm.html) 
  + [`PoissonFlowFromFlux`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromflux.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxDCTI`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxdcti.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxType1`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxtype1.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxWithAxes`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/poissonflowfromfluxwithaxes.html) We will treat the first dimension as `x' and the second
  + [`addForcingFluxesForReservoirGroupAtTime`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/addforcingfluxesforreservoirgroupattime.html) 
  + [`areEnergyReservoirsComplete`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/areenergyreservoirscomplete.html) 
  + [`areTriadComponentsComplete`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/aretriadcomponentscomplete.html) 
  + [`diagfile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagfile.html) 
  + [`diagnosticsHasExplicitAntialiasing`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagnosticshasexplicitantialiasing.html) 
  + [`diagpath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/diagpath.html) 
  + [`escale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/escale.html) 
  + [`escale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/escale_units.html) 
  + [`fancyNameForName`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/fancynameforname.html) 
  + [`filterFluxesForReservoir`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/filterfluxesforreservoir.html) 
  + [`flux_scale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/flux_scale.html) 
  + [`flux_scale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/flux_scale_units.html) 
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
  + [`logWavelengthAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/logwavelengthaxis.html) To use this:
  + [`omegaAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/omegaaxis.html) 
  + [`omega_jk`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/omega_jk.html) 
  + [`overlayFrequencyContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlayfrequencycontours.html) 
  + [`overlayGeostrophicKineticPotentialFractionContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlaygeostrophickineticpotentialfractioncontours.html) 
  + [`overlayGeostrophicKineticPotentialRatioContours`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/overlaygeostrophickineticpotentialratiocontours.html) 
  + [`plotPoissonFlowOverPcolor`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/plotpoissonflowoverpcolor.html) 
  + [`pseudoRadialBinning`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/pseudoradialbinning.html) 
  + [`setLogWavelengthXAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/setlogwavelengthxaxis.html) 
  + [`showRossbyRadiusYAxis`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/showrossbyradiusyaxis.html) 
  + [`symmetricTintMap`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/symmetrictintmap.html) Colormap going from color -> tinted white -> color
  + [`t_diag`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/t_diag.html) 
  + [`t_wv`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/t_wv.html) 
  + [`tscale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/tscale.html) Default scaling and units for time, energy, and flux
  + [`tscale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/tscale_units.html) 
  + [`version`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/version.html) Locate folder containing this class file
  + [`wvaapath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvaapath.html) 
  + [`wvfile`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvfile.html) 
  + [`wvpath`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvpath.html) 
  + [`wvt`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvt.html) 
  + [`wvt_aa`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/wvt_aa.html) 
  + [`z_flux_scale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/z_flux_scale.html) 
  + [`z_flux_scale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/z_flux_scale_units.html) 
  + [`zscale`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/zscale.html) 
  + [`zscale_units`](/wave-vortex-model-diagnostics/classes/wvdiagnostics/zscale_units.html) 


---