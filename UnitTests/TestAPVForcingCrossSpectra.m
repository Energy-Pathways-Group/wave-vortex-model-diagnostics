% basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = "/Volumes/Samsung_T7/CimRuns_June2025/output/";
% basedir = '/Users/cwortham/Documents/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns/output/';
% basedir = '/Volumes/SanDiskExtremePro/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns_June2025_v2/output/';

runNumber=9; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");
% wvd = WVDiagnostics(basedir + getRunParameters(runNumber) + ".nc");

% wvt = wvd.wvt;
% ncfile = wvd.wvfile;

%%


%%

options.shouldMeasureAntialiasingFlux = true;

if options.shouldMeasureAntialiasingFlux
    [wvt_lowres, ncfile] = WVTransform.waveVortexTransformFromFile(wvd.wvfile.path,iTime=Inf);
    wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
else
    [wvt, ncfile] = WVTransform.waveVortexTransformFromFile(wvd.wvfile.path,iTime=Inf);
end

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APEOperation(wvt));
wvt.addOperation(APVOperation());
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

forcingNames = wvt.forcingNames;

%%
indices = 51;
Z2_qgpv_t = zeros(length(indices),1);
Z2_t = zeros(length(indices),1);
iIndex = 1;

if options.shouldMeasureAntialiasingFlux
    wvt_lowres.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
    wvt.t = wvt_lowres.t;
    [wvt.A0,wvt.Ap,wvt.Am] = wvt_lowres.spectralVariableWithResolution(wvt,wvt_lowres.A0,wvt_lowres.Ap,wvt_lowres.Am);
else
    wvt.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
end

% wvt.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
F = wvt.fluxForForcing();

fprintf("Enstrophy forcing:\n");
enstrophyScale = wvt.f*wvt.f/(86400*365);
eta_true = wvt.eta_true;
for iForce=1:1%length(forcingNames)
    if isa(wvt,"WVTransformHydrostatic")
        [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        if forcingNames(iForce) == "nonlinear advection"
            Fu = Fu + wvt.f*wvt.v;
            Fv = Fv - wvt.f*wvt.u;
        end
        DF_x = - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    else
        [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        if forcingNames(iForce) == "nonlinear advection"
            Fu = Fu + wvt.f*wvt.v;
            Fv = Fv - wvt.f*wvt.u;
        end
        DF_x = wvt.diffY(Fw) - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu) - wvt.diffX(Fw);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    end
    Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(iForce)}.F0 .* conj(wvt.A0) );
    if isscalar(indices)
        fprintf("Total quadratic enstrophy for " + forcingNames(iForce) + " forcing: " + sum(Z(:))/enstrophyScale + "\n");
    end
    Z2_qgpv_t(iIndex) = sum(Z(:))/enstrophyScale;

    % F_pv = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(Feta);
    % Z2 = wvt.qgpv .* F_pv;
    % if isscalar(indices)
    %     fprintf("Total approximate enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
    % end
    % F0_pv = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(F_pv));


    % zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
    % zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
    % zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

    if forcingNames(iForce) == "nonlinear advection"
        G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*(Feta + wvt.w);
        % G_eta = (-wvt.u .* wvt.diffX(eta_true) - wvt.v .* wvt.diffY(eta_true) - wvt.w .* (wvt.diffZG(eta_true) - 1));
    else
        G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*Feta;
    end
    % G_eta = Feta;
    FZ_L = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(G_eta);
    FZ_NL = - wvt.zeta_x .* wvt.diffX(G_eta) - wvt.zeta_y .* wvt.diffY(G_eta)- wvt.zeta_z .* wvt.diffZG(G_eta);
    FZ_NL = FZ_NL - wvt.diffX(eta_true) .* DF_x - wvt.diffY(eta_true) .* DF_y - wvt.diffZG(eta_true) .* DF_z;
    FZ = FZ_L + FZ_NL;
    Z2 = wvt.apv .* FZ;

    Z2_t(iIndex) = int_vol(Z2)/enstrophyScale;

    S_f = wvd.crossSpectrumWithFgTransform(FZ,wvt.apv);

    if isscalar(indices)
        fprintf("Total nonlinear enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + " (" + sum(S_f(:))/enstrophyScale + ")\n");
    end
    
end