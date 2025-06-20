basedir = "/Users/Shared/CimRuns_June2025/output/";
runNumber = 1;
additionalDays = 1000;

runName = getRunParameters(runNumber);

filepath = basedir + runName + ".nc";

wvt = WVTransform.waveVortexTransformFromFile(filepath);

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APEOperation(wvt));
wvt.addOperation(APVOperation());
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

F = wvt.fluxForForcing();
forcingNames = wvt.forcingNames;

%%
fprintf("Enstrophy forcing:\n");
enstrophyScale = wvt.f*wvt.f/(86400*365);
for iForce=1:length(forcingNames)
    [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));

    Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(iForce)}.F0 .* conj(wvt.A0) );
    fprintf("Total quadratic enstrophy for " + forcingNames(iForce) + " forcing: " + sum(Z(:))/enstrophyScale + "\n");

    Z2 = wvt.qgpv .* (wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(Feta));
    fprintf("Total approximate enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
end