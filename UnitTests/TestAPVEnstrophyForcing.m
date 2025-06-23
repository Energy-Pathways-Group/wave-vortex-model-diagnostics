basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
runNumber = 1;

runName = getRunParameters(runNumber);

filepath = basedir + runName + ".nc";

[wvt, ncfile] = WVTransform.waveVortexTransformFromFile(filepath);

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APEOperation(wvt));
wvt.addOperation(APVOperation());
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

forcingNames = wvt.forcingNames;

%%

wvt.initFromNetCDFFile(ncfile,iTime=870);
F = wvt.fluxForForcing();

%%

fprintf("Enstrophy forcing:\n");
enstrophyScale = wvt.f*wvt.f/(86400*365);
eta_true = wvt.eta_true;
for iForce=1:length(forcingNames)
    if false %isa(wvt,"WVTransformHydrostatic")
        [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        DF_x = - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    else

        [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        DF_x = wvt.diffY(Fw) - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu) - wvt.diffX(Fw);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    end
    Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(iForce)}.F0 .* conj(wvt.A0) );
    fprintf("Total quadratic enstrophy for " + forcingNames(iForce) + " forcing: " + sum(Z(:))/enstrophyScale + "\n");

    F_pv = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(Feta);
    Z2 = wvt.qgpv .* F_pv;
    fprintf("Total approximate enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
    F0_pv = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(F_pv));


    % zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
    % zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
    % zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

    G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*Feta;
    % G_eta = Feta;
    Z_NL = - wvt.zeta_x .* wvt.diffX(G_eta) - wvt.zeta_y .* wvt.diffY(G_eta)- wvt.zeta_z .* wvt.diffZG(G_eta);
    Z_NL = Z_NL - wvt.diffX(eta_true) .* DF_x - wvt.diffY(eta_true) .* DF_y - wvt.diffZG(eta_true) .* DF_z;
    Z2 = wvt.apv .* (wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(G_eta) + Z_NL);

    fprintf("Total nonlinear enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
end