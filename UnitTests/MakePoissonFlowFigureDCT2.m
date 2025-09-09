basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=9; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
energy_fluxes = wvd.exactEnergyFluxesTemporalAverage(timeIndices=51:251);
enstrophy_fluxes = wvd.exactEnstrophyFluxesTemporalAverage(timeIndices=51:251);

inertial_fluxes = wvd.quadraticEnergyTriadFluxesTemporalAverage(timeIndices=51:251);
ggg = inertial_fluxes(1).te_gmda;
wgg = inertial_fluxes(1).te_wave + inertial_fluxes(2).te_gmda + inertial_fluxes(3).te_gmda;
wwg = inertial_fluxes(2).te_wave + inertial_fluxes(3).te_wave + inertial_fluxes(4).te_gmda;
www = inertial_fluxes(4).te_wave;

%%
flux = energy_fluxes(1).te/wvd.flux_scale;
flux = enstrophy_fluxes(1).Z0/wvd.z_flux_scale;
flux = www;

[X,Y,U,V] = wvd.PoissonFlowFromFlux(flux.');
% figure
% quiver(X,Y,U,V,Color=0*[1 1 1])


kRadial = wvd.kRadial;
% kRadial = kRadial + (kRadial(2)-kRadial(1))/2;
% kRadial(1) = 0.5*wvd.kRadial(2);
jWavenumber = wvd.jWavenumber;
% jWavenumber = jWavenumber + (jWavenumber(2)-jWavenumber(1))/2;

% jWavenumber(1) = 0.5*wvd.jWavenumber(2);

figure, pcolor(kRadial,jWavenumber,flux); shading flat;

% pbaspect([wvd.jWavenumber(2)/wvd.kRadial(2) 1 1])

% set(gca,'YDir','reverse')
% set(gca,'XDir','reverse')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
xlim([kRadial(1) 1.6e-3])
ylim([jWavenumber(1) 1.6e-3])

colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(flux(:)))*[-1 1]/2)
colorbar("eastoutside")
hold on,
quiver(X,Y,U,V,Color=0*[1 1 1])

%%
figure
plot(wvd.kRadial,cumsum(squeeze(sum(flux,1)))), hold on
plot(wvd.kRadial,-(squeeze(sum(U/wvd.kRadial(2),2))))

%%
figure
plot(wvd.jWavenumber,cumsum(squeeze(sum(flux,2)))), hold on
plot(wvd.jWavenumber,-(squeeze(sum(U/wvd.jWavenumber(2),1))))

%%

[X,Y,U,V] = wvd.PoissonFlowFromFluxType1(flux.');
figure
quiver(X,Y,U,V,Color=0*[1 1 1])