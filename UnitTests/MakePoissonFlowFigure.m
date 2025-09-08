% basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=18; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
energy_fluxes = wvd.exactEnergyFluxesTemporalAverage();


%%
flux = energy_fluxes(1).te;
jWavenumber = 1./sqrt(wvd.Lr2);
jWavenumber(1) = 0; % barotropic mode is a mean?
[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvd.kRadial,jWavenumber,flux.');

figure, jpcolor(wvd.kRadial,jWavenumber,flux); shading flat;
colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(flux(:)))*[-1 1])
colorbar("eastoutside")
hold on,
quiver(X,Y,10*U,10*V,Color=0*[1 1 1])
% xlim([0 5e-4])
% ylim([0 5e-4])