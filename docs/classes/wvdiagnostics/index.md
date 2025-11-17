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
      + [`CosineTransformBack`](/classes/wvdiagnostics/cosinetransformback.html) Fast Discrete Inverse Cosine Transform.
      + [`CosineTransformForward`](/classes/wvdiagnostics/cosinetransformforward.html) Fast Discrete Cosine Transform (DCT-I).
      + [`SineTransformBack`](/classes/wvdiagnostics/sinetransformback.html) Fast Discrete Inverse Sine Transform.
      + [`transformToKePeAxis`](/classes/wvdiagnostics/transformtokepeaxis.html) Transform To Ke Pe Axis.
      + [`transformToOmegaAxis`](/classes/wvdiagnostics/transformtoomegaaxis.html) transforms in the from (j,kRadial) to omegaAxis.
      + [`transformToPseudoRadialWavenumber`](/classes/wvdiagnostics/transformtopseudoradialwavenumber.html) transforms in the from (j,kRadial) to kPseudoRadial.
      + [`transformToPseudoRadialWavenumberA0`](/classes/wvdiagnostics/transformtopseudoradialwavenumbera0.html) transforms in the from (j,kRadial) to kPseudoRadial.
      + [`transformToPseudoRadialWavenumberApm`](/classes/wvdiagnostics/transformtopseudoradialwavenumberapm.html) transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial.
  + Spectral
    + General
      + [`SineTransformForward`](/classes/wvdiagnostics/sinetransformforward.html) Fast Discrete Sine Transform (DST-I).
+ Diagnostics
  + General — Misc
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`DCT1`](/classes/wvdiagnostics/dct1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`DCT2`](/classes/wvdiagnostics/dct2.html) CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix.
      + [`DST1`](/classes/wvdiagnostics/dst1.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`DST2`](/classes/wvdiagnostics/dst2.html) SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix.
      + [`cmocean`](/classes/wvdiagnostics/cmocean.html) returns perceptually-uniform colormaps created by Kristen Thyng.
      + [`geostrophicGeostrophicWaveEnergy`](/classes/wvdiagnostics/geostrophicgeostrophicwaveenergy.html) Note that.
      + [`iDCT1`](/classes/wvdiagnostics/idct1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`iDCT2`](/classes/wvdiagnostics/idct2.html) InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix.
      + [`iDST1`](/classes/wvdiagnostics/idst1.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
      + [`iDST2`](/classes/wvdiagnostics/idst2.html) InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix.
      + [`waveWaveGeostrophicEnergy`](/classes/wvdiagnostics/wavewavegeostrophicenergy.html) Note that.
      + [`waveWaveGeostrophicEnergyForMode`](/classes/wvdiagnostics/wavewavegeostrophicenergyformode.html) Note that.
    + Fluxes over time, [t 1]
      + [`quadraticEnergyMirrorTriadsUndamped`](/classes/wvdiagnostics/quadraticenergymirrortriadsundamped.html) Quadratic Energy Mirror Triads Undamped.
  + Flux diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`RescalePoissonFlowFluxForLogSpace`](/classes/wvdiagnostics/rescalepoissonflowfluxforlogspace.html) Rescale Poisson Flow Flux For Log Space.
  + Flux diagnostics
    + Triad interactions
      + [`addTriadFluxesForReservoirGroupAtTime`](/classes/wvdiagnostics/addtriadfluxesforreservoirgroupattime.html) Add Triad Fluxes For Reservoir Group At Time.
    + General
      + [`create1DMirrorFluxes`](/classes/wvdiagnostics/create1dmirrorfluxes.html) Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
      + [`create2DMirrorFluxes`](/classes/wvdiagnostics/create2dmirrorfluxes.html) Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
      + [`fluxesForReservoirGroup`](/classes/wvdiagnostics/fluxesforreservoirgroup.html) Fluxes For Reservoir Group.
      + [`showDampingFluxVsPseudolength`](/classes/wvdiagnostics/showdampingfluxvspseudolength.html) Show Damping Flux Vs Pseudolength.
  + General
    + Misc
      + [`createDiagnosticsFile`](/classes/wvdiagnostics/creatediagnosticsfile.html) Create a new diagnostics file and compute diagnostics from WVModel output.
  + Enstrophy fluxes
    + General
      + [`createEnstrophyFluxSummaryTable`](/classes/wvdiagnostics/createenstrophyfluxsummarytable.html) Create a LaTeX table summarizing enstrophy source/sink fluxes.
  + Spectra — Cross-spectra
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`crossSpectrumWithFgTransform`](/classes/wvdiagnostics/crossspectrumwithfgtransform.html) Cross Spectrum With Fg Transform.
      + [`crossSpectrumWithGgTransform`](/classes/wvdiagnostics/crossspectrumwithggtransform.html) Cross Spectrum With Gg Transform.
  + Energy fluxes — General
    + Fluxes, [j kRadial t]
      + [`exactEnergyFluxes`](/classes/wvdiagnostics/exactenergyfluxes.html) Return the exact energy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnergyFluxes`](/classes/wvdiagnostics/quadraticenergyfluxes.html) Return the energy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`exactEnergyFluxesOverTime`](/classes/wvdiagnostics/exactenergyfluxesovertime.html) Exact energy fluxes over time.
      + [`quadraticEnergyFluxesOverTime`](/classes/wvdiagnostics/quadraticenergyfluxesovertime.html) Compute forcing fluxes over time.
  + Energy fluxes — Time/space averages
    + Flux averages, scalar [1 1]
      + [`exactEnergyFluxesSpatialTemporalAverage`](/classes/wvdiagnostics/exactenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact forcing fluxes.
      + [`quadraticEnergyFluxesSpatialTemporalAverage`](/classes/wvdiagnostics/quadraticenergyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of forcing fluxes.
    + Fluxes in space, [j kRadial]
      + [`exactEnergyFluxesTemporalAverage`](/classes/wvdiagnostics/exactenergyfluxestemporalaverage.html) Temporally averaged exact energy fluxes.
      + [`quadraticEnergyFluxesTemporalAverage`](/classes/wvdiagnostics/quadraticenergyfluxestemporalaverage.html) Compute temporally averaged forcing fluxes.
  + Energy
    + Time series
      + [`exactEnergyOverTime`](/classes/wvdiagnostics/exactenergyovertime.html) Exact Energy Over Time.
      + [`quadraticEnergyOverTime`](/classes/wvdiagnostics/quadraticenergyovertime.html) Compute energy for each reservoir over time.
  + Enstrophy fluxes — General
    + Fluxes, [j kRadial t]
      + [`exactEnstrophyFluxes`](/classes/wvdiagnostics/exactenstrophyfluxes.html) Return the available potential enstrophy flux from the forcing terms, [j kRadial t].
      + [`quadraticEnstrophyFluxes`](/classes/wvdiagnostics/quadraticenstrophyfluxes.html) Return the enstrophy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`exactEnstrophyFluxesOverTime`](/classes/wvdiagnostics/exactenstrophyfluxesovertime.html) Compute exact enstrophy fluxes over time.
      + [`quadraticEnstrophyFluxesOverTime`](/classes/wvdiagnostics/quadraticenstrophyfluxesovertime.html) Compute enstrophy fluxes over time.
  + Enstrophy fluxes — Time/space averages
    + Flux averages, scalar [1 1]
      + [`exactEnstrophyFluxesSpatialTemporalAverage`](/classes/wvdiagnostics/exactenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the exact enstrophy fluxes.
      + [`quadraticEnstrophyFluxesSpatialTemporalAverage`](/classes/wvdiagnostics/quadraticenstrophyfluxesspatialtemporalaverage.html) Compute spatial-temporal average of the qgpv enstrophy fluxes.
    + Fluxes in space, [j kRadial]
      + [`exactEnstrophyFluxesTemporalAverage`](/classes/wvdiagnostics/exactenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
      + [`quadraticEnstrophyFluxesTemporalAverage`](/classes/wvdiagnostics/quadraticenstrophyfluxestemporalaverage.html) Compute temporally averaged enstrophy fluxes.
  + Enstrophy
    + Time series
      + [`exactEnstrophyOverTime`](/classes/wvdiagnostics/exactenstrophyovertime.html) Exact Enstrophy Over Time.
      + [`quadraticEnstrophyOverTime`](/classes/wvdiagnostics/quadraticenstrophyovertime.html) Quadratic Enstrophy Over Time.
  + Energy fluxes — Triad interactions — Mirror pairs
    + Fluxes in space, [sparseKRadialAxis 1]
      + [`quadraticEnergyMirrorTriadFluxes1D`](/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d.html) Quadratic Energy Mirror Triad Fluxes1 D.
      + [`quadraticEnergyMirrorTriadFluxes1D_kepe`](/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_kepe.html) Quadratic Energy Mirror Triad Fluxes1 D kepe.
      + [`quadraticEnergyMirrorTriadFluxes1D_omega`](/classes/wvdiagnostics/quadraticenergymirrortriadfluxes1d_omega.html) Quadratic Energy Mirror Triad Fluxes1 D omega.
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyMirrorTriadFluxes2D`](/classes/wvdiagnostics/quadraticenergymirrortriadfluxes2d.html) Quadratic Energy Mirror Triad Fluxes2 D.
  + Energy fluxes — Triad interactions — Primary
    + Fluxes in space, [sparseKRadialAxis 1]
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D`](/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe`](/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_kepe.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D kepe.
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega`](/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage1d_omega.html) Quadratic Energy Primary Triad Fluxes Temporal Average1 D omega.
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyPrimaryTriadFluxesTemporalAverage2D`](/classes/wvdiagnostics/quadraticenergyprimarytriadfluxestemporalaverage2d.html) outputGrid determines whether or not the fluxes get downsampled to the.
  + Energy fluxes — Triad interactions
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`quadraticEnergyTriadFluxWWGWave`](/classes/wvdiagnostics/quadraticenergytriadfluxwwgwave.html) Quadratic Energy Triad Flux WWGWave.
    + Fluxes, [j kRadial t]
      + [`quadraticEnergyTriadFluxes`](/classes/wvdiagnostics/quadraticenergytriadfluxes.html) Return the energy flux from the inertial terms, specified as triad components.
    + Fluxes over time, [t 1]
      + [`quadraticEnergyTriadFluxesOverTime`](/classes/wvdiagnostics/quadraticenergytriadfluxesovertime.html) Compute inertial fluxes over time.
    + Flux averages, scalar [1 1]
      + [`quadraticEnergyTriadFluxesSpatialTemporalAverage`](/classes/wvdiagnostics/quadraticenergytriadfluxesspatialtemporalaverage.html) Compute spatial-temporal average of inertial fluxes.
    + Fluxes in space, [j kRadial]
      + [`quadraticEnergyTriadFluxesTemporalAverage`](/classes/wvdiagnostics/quadraticenergytriadfluxestemporalaverage.html) Computes the temporally averaged quadratic triad fluxes.
  + Enstrophy fluxes — Triad interactions
    + Fluxes, [j kRadial t]
      + [`quadraticEnstrophyTriadFluxes`](/classes/wvdiagnostics/quadraticenstrophytriadfluxes.html) Return the enstrophy flux from the forcing terms.
    + Fluxes over time, [t 1]
      + [`quadraticEnstrophyTriadFluxesOverTime`](/classes/wvdiagnostics/quadraticenstrophytriadfluxesovertime.html) Compute enstrophy inertial (aka, triad) fluxes over time.
  + Spectra — Potential energy
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`spectrumWithFgTransform`](/classes/wvdiagnostics/spectrumwithfgtransform.html) Spectrum With Fg Transform.
      + [`spectrumWithGgTransform`](/classes/wvdiagnostics/spectrumwithggtransform.html) Spectrum With Gg Transform.
+ Constructor
  + [`WVDiagnostics`](/classes/wvdiagnostics/wvdiagnostics.html) Initializes the WVDiagnostics object, loads the wave-vortex transform and diagnostics files.
+ Utilities
  + Colormaps — Crameri
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`crameri`](/classes/wvdiagnostics/crameri.html) returns perceptually-uniform scientific colormaps created.
  + Sparse matrices — Axis binning
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`sparseJKAxisBinMatrices`](/classes/wvdiagnostics/sparsejkaxisbinmatrices.html) Sparse JKAxis Bin Matrices.
      + [`sparseJWavenumberAxis`](/classes/wvdiagnostics/sparsejwavenumberaxis.html) Sparse JWavenumber Axis.
      + [`sparseKRadialAxis`](/classes/wvdiagnostics/sparsekradialaxis.html) Sparse KRadial Axis.
      + [`sparseKePeAxis`](/classes/wvdiagnostics/sparsekepeaxis.html) Sparse Ke Pe Axis.
      + [`sparseOmegaAxis`](/classes/wvdiagnostics/sparseomegaaxis.html) Sparse Omega Axis.
      + [`sparsePseudoRadialAxis`](/classes/wvdiagnostics/sparsepseudoradialaxis.html) Sparse Pseudo Radial Axis.
+ Configuration
  + [`setEnergyUnits`](/classes/wvdiagnostics/setenergyunits.html) Set the time and energy scaling and units for plotting and output.
  + Reservoirs
    + Grouping
      + [`createReservoirGroup`](/classes/wvdiagnostics/createreservoirgroup.html) Create a reservoir group in the diagnostics NetCDF and populate it.
      + [`filterEnergyForSourcesSinksReservoirs`](/classes/wvdiagnostics/filterenergyforsourcessinksreservoirs.html) This function returns values assuming three reservoirs: geo, wave, and.
      + [`filterEnergyForSourcesSinksReservoirsOld`](/classes/wvdiagnostics/filterenergyforsourcessinksreservoirsold.html) This function returns values assuming three reservoirs: geo, wave, and.
      + [`summarizeSourcesSinksReservoirs`](/classes/wvdiagnostics/summarizesourcessinksreservoirs.html) Summarize Sources Sinks Reservoirs.
      + [`variablesForReservoirGroup`](/classes/wvdiagnostics/variablesforreservoirgroup.html) Variables For Reservoir Group.
+ Figures
  + Time series
    + Diagnostics
      + [`plotEnergyFluxOverTime`](/classes/wvdiagnostics/plotenergyfluxovertime.html) Plot energy fluxes for reservoirs as a function of time.
      + [`plotEnergyOverTime`](/classes/wvdiagnostics/plotenergyovertime.html) Plot energy for each reservoir over time.
      + [`plotEnergyTriadFluxOverTime`](/classes/wvdiagnostics/plotenergytriadfluxovertime.html) Plot inertial flux for each reservoir over time.
      + [`plotEnstrophyFluxOverTime`](/classes/wvdiagnostics/plotenstrophyfluxovertime.html) Plot enstrophy fluxes for reservoirs as a function of time.
      + [`plotEnstrophyOverTime`](/classes/wvdiagnostics/plotenstrophyovertime.html) Plot enstrophy over time.
      + [`plotEnstrophyTriadFluxOverTime`](/classes/wvdiagnostics/plotenstrophytriadfluxovertime.html) Plot enstrophy triad fluxes over time.
  + Diagnostics
    + General
      + [`plotEnergyFluxTemporalAverage`](/classes/wvdiagnostics/plotenergyfluxtemporalaverage.html) Plot temporally averaged energy flux diagnostics.
      + [`plotEnergyFluxes1D`](/classes/wvdiagnostics/plotenergyfluxes1d.html) Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency.
      + [`plotEnstrophyFluxTemporalAverage`](/classes/wvdiagnostics/plotenstrophyfluxtemporalaverage.html) Plot temporally averaged enstrophy flux diagnostics.
      + [`plotFluidDecompositionMultipanel`](/classes/wvdiagnostics/plotfluiddecompositionmultipanel.html) Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
      + [`plotFluidStateMultipanel`](/classes/wvdiagnostics/plotfluidstatemultipanel.html) Plot multipanel summary of fluid state and spectra
      + [`plotSourcesSinksForReservoirGroup`](/classes/wvdiagnostics/plotsourcessinksforreservoirgroup.html) Plot sources, sinks, and reservoirs diagram for a reservoir group.
      + [`plotSourcesSinksReservoirsDiagram`](/classes/wvdiagnostics/plotsourcessinksreservoirsdiagram.html) Plot sources, sinks, and reservoirs diagram.
      + [`plotSourcesSinksReservoirsDiagramWithClosureRegion`](/classes/wvdiagnostics/plotsourcessinksreservoirsdiagramwithclosureregion.html) Plot sources, sinks, and reservoirs diagram.
  + Spectra
    + Potential energy
      + [`plotEnergySpectrum`](/classes/wvdiagnostics/plotenergyspectrum.html) Plot the wave/geostrophic energy spectra at a given time.
      + [`plotEnergySpectrumOLD`](/classes/wvdiagnostics/plotenergyspectrumold.html) Plot the wave/geostrophic energy spectra at a given time.
      + [`plotEnstrophySpectrum`](/classes/wvdiagnostics/plotenstrophyspectrum.html) Plot the enstrophy spectrum at a given time.
      + [`plotEnstrophySpectrumOLD`](/classes/wvdiagnostics/plotenstrophyspectrumold.html) Plot the geostrophic enstrophy spectrum at a given time.
      + [`plotMooringRotarySpectrum`](/classes/wvdiagnostics/plotmooringrotaryspectrum.html) Plot rotary spectra from mooring velocity time series.
      + [`plotPotentialEnergySpectrum`](/classes/wvdiagnostics/plotpotentialenergyspectrum.html) Plot energy for each reservoir over time.
  + Diagnostics — General
    + Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
      + [`plotPoissonFlowOverContours`](/classes/wvdiagnostics/plotpoissonflowovercontours.html) Plot Poisson Flow Over Contours.
+ Other
  + [`Lr2`](/classes/wvdiagnostics/lr2.html) 
  + [`Lr2_pm`](/classes/wvdiagnostics/lr2_pm.html) 
  + [`PoissonFlowFromFlux`](/classes/wvdiagnostics/poissonflowfromflux.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxDCTI`](/classes/wvdiagnostics/poissonflowfromfluxdcti.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxType1`](/classes/wvdiagnostics/poissonflowfromfluxtype1.html) We will treat the first dimension as `x' and the second
  + [`PoissonFlowFromFluxWithAxes`](/classes/wvdiagnostics/poissonflowfromfluxwithaxes.html) We will treat the first dimension as `x' and the second
  + [`addForcingFluxesForReservoirGroupAtTime`](/classes/wvdiagnostics/addforcingfluxesforreservoirgroupattime.html) 
  + [`areEnergyReservoirsComplete`](/classes/wvdiagnostics/areenergyreservoirscomplete.html) 
  + [`areTriadComponentsComplete`](/classes/wvdiagnostics/aretriadcomponentscomplete.html) 
  + [`diagfile`](/classes/wvdiagnostics/diagfile.html) 
  + [`diagnosticsHasExplicitAntialiasing`](/classes/wvdiagnostics/diagnosticshasexplicitantialiasing.html) 
  + [`diagpath`](/classes/wvdiagnostics/diagpath.html) 
  + [`escale`](/classes/wvdiagnostics/escale.html) 
  + [`escale_units`](/classes/wvdiagnostics/escale_units.html) 
  + [`fancyNameForName`](/classes/wvdiagnostics/fancynameforname.html) 
  + [`filterFluxesForReservoir`](/classes/wvdiagnostics/filterfluxesforreservoir.html) 
  + [`flux_scale`](/classes/wvdiagnostics/flux_scale.html) 
  + [`flux_scale_units`](/classes/wvdiagnostics/flux_scale_units.html) 
  + [`forcingNames`](/classes/wvdiagnostics/forcingnames.html) 
  + [`geo_hke_jk`](/classes/wvdiagnostics/geo_hke_jk.html) 
  + [`geo_pe_jk`](/classes/wvdiagnostics/geo_pe_jk.html) 
  + [`iTime`](/classes/wvdiagnostics/itime.html) 
  + [`iTimeChanged`](/classes/wvdiagnostics/itimechanged.html) 
  + [`j`](/classes/wvdiagnostics/j.html) 
  + [`jWavenumber`](/classes/wvdiagnostics/jwavenumber.html) 
  + [`kPseudoRadial`](/classes/wvdiagnostics/kpseudoradial.html) 
  + [`kRadial`](/classes/wvdiagnostics/kradial.html) 
  + [`kePeAxis`](/classes/wvdiagnostics/kepeaxis.html) 
  + [`logWavelengthAxis`](/classes/wvdiagnostics/logwavelengthaxis.html) To use this:
  + [`omegaAxis`](/classes/wvdiagnostics/omegaaxis.html) 
  + [`omega_jk`](/classes/wvdiagnostics/omega_jk.html) 
  + [`overlayFrequencyContours`](/classes/wvdiagnostics/overlayfrequencycontours.html) 
  + [`overlayGeostrophicKineticPotentialFractionContours`](/classes/wvdiagnostics/overlaygeostrophickineticpotentialfractioncontours.html) 
  + [`overlayGeostrophicKineticPotentialRatioContours`](/classes/wvdiagnostics/overlaygeostrophickineticpotentialratiocontours.html) 
  + [`plotPoissonFlowOverPcolor`](/classes/wvdiagnostics/plotpoissonflowoverpcolor.html) 
  + [`pseudoRadialBinning`](/classes/wvdiagnostics/pseudoradialbinning.html) 
  + [`setLogWavelengthXAxis`](/classes/wvdiagnostics/setlogwavelengthxaxis.html) 
  + [`showRossbyRadiusYAxis`](/classes/wvdiagnostics/showrossbyradiusyaxis.html) 
  + [`symmetricTintMap`](/classes/wvdiagnostics/symmetrictintmap.html) Colormap going from color -> tinted white -> color
  + [`t_diag`](/classes/wvdiagnostics/t_diag.html) 
  + [`t_wv`](/classes/wvdiagnostics/t_wv.html) 
  + [`tscale`](/classes/wvdiagnostics/tscale.html) Default scaling and units for time, energy, and flux
  + [`tscale_units`](/classes/wvdiagnostics/tscale_units.html) 
  + [`version`](/classes/wvdiagnostics/version.html) Locate folder containing this class file
  + [`wvaapath`](/classes/wvdiagnostics/wvaapath.html) 
  + [`wvfile`](/classes/wvdiagnostics/wvfile.html) 
  + [`wvpath`](/classes/wvdiagnostics/wvpath.html) 
  + [`wvt`](/classes/wvdiagnostics/wvt.html) 
  + [`wvt_aa`](/classes/wvdiagnostics/wvt_aa.html) 
  + [`z_flux_scale`](/classes/wvdiagnostics/z_flux_scale.html) 
  + [`z_flux_scale_units`](/classes/wvdiagnostics/z_flux_scale_units.html) 
  + [`zscale`](/classes/wvdiagnostics/zscale.html) 
  + [`zscale_units`](/classes/wvdiagnostics/zscale_units.html) 


---