basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = "/Users/Shared/CimRuns_June2025/output/";
wvd = WVDiagnostics(basedir + replace(getRunParameters(9),"256","256") + ".nc");
wvt = wvd.wvt;
timeIndices = 2751:3001;

%%
timeIndices = 2951:3001;
wvd.createReservoirGroup(timeIndices=2951:3001);

%%
[transferFlux, forcingFlux] = wvd.fluxesForReservoirGroup(timeIndices=timeIndices);
