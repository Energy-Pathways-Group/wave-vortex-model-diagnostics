basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

%%
runNumber=1; runName = "hydrostatic: geostrophic";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
runNumber=9; runName = "hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvd.plotFluidStateMultipanel();

%%
wvd.plotFluidDecompositionMultipanel();

%%
wvd.plotEnstrophySpectrum();

%%
wvd.plotEnergySpectrum();

%%
wvd.plotMooringRotarySpectrum();

%%
wvd.plotEnergyFluxTemporalAverage(approximation="exact",overSaturationFactor=100);

%%
wvd.plotEnergyFluxTemporalAverage(approximation="quadratic");

%%
wvd.plotEnstrophyOverTime();

%%
wvd.plotEnergyOverTime();

%%
wvd.plotEnergyFluxOverTime(approximation="exact",filter=@(v,t) movmean(v,21));

%%
wvd.plotEnergyFluxOverTime(approximation="quadratic",filter=@(v,t) movmean(v,21));

%%
wvd.plotEnergyTriadFluxOverTime(filter=@(v,t) movmean(v,21));

%%
wvd.plotEnstrophyFluxOverTime(filter=@(v,t) movmean(v,21));

%%
wvd.plotExactEnstrophyFluxOverTime(filter=@(v,t) movmean(v,21));