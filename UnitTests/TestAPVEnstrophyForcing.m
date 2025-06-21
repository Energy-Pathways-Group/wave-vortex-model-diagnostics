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

    zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
    zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
    zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

    % DF_x = wvt.diffY(Fw) - wvt.diffZF(Fv); % w_y - v_z
    % DF_y = wvt.diffZF(Fu) - wvt.diffX(Fw;  % u_z - w_x
    % DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y

    DF_x = - wvt.diffZF(Fv); % w_y - v_z
    DF_y = wvt.diffZF(Fu);  % u_z - w_x
    DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y

    Z_NL = - zeta_x .* wvt.diffX(Feta) - zeta_y .* wvt.diffY(Feta)- zeta_z .* wvt.diffZG(Feta);
    Z_NL = Z_NL - wvt.diffX(eta_true) .* DF_x - wvt.diffY(eta_true) .* DF_y - wvt.diffZG(eta_true) .* DF_z; 
end