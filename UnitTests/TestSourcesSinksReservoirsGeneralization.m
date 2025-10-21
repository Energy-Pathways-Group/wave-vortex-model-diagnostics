% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
basedir = "/Users/Shared/CimRuns_June2025/output/";
wvd = WVDiagnostics(basedir + replace(getRunParameters(18),"256","512") + ".nc");
wvt = wvd.wvt;
% timeIndices = 2751:3001;
timeIndices = 51:251;

%%
% timeIndices = 51:251;
% wvd.createReservoirGroup(timeIndices=timeIndices);

%%
[transferFlux, forcingFlux, ddt, energy] = wvd.fluxesForReservoirGroup(timeIndices=51:251);
