basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

%%
runNumber=1; runName = "hydrostatic: geostrophic";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
runNumber=9; runName = "hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvt = wvd.wvt;

U2 = wvd.spectrumWithFgTransform(wvt.u);
V2 = wvd.spectrumWithFgTransform(wvt.v);
N2 = wvd.spectrumWithGgTransform(wvt.eta);
% N2 = wvd.crossSpectrumWithGgTransform(wvt.eta,wvt.eta);

E = 0.5*(U2 + V2 + N2);

E2 = wvt.Apm_TE_factor.*( abs(wvt.Ap).^2 + abs(wvt.Am).^2 )  + wvt.A0_TE_factor.*( abs(wvt.A0).^2);

max(abs(E(:)-E2(:)))

%%
runNumber=18; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvt = wvd.wvt;

U2 = wvd.spectrumWithFgTransform(wvt.u);
V2 = wvd.spectrumWithFgTransform(wvt.v);
W2 = wvd.spectrumWithGgTransform(wvt.w);
rho_nm = @(z) wvt.rhoFunction(z) - wvt.rho0;
rho_nm_z = @(z) -(wvt.rho0/wvt.g)*wvt.N2Function(z);
eta_nl = (wvt.g/wvt.rho0)*rho_nm(wvt.Z - wvt.eta_true) ./ wvt.N2Function(wvt.Z);
N2 = wvd.crossSpectrumWithGgTransform(wvt.eta_true,eta_nl);

E = 0.5*(U2 + V2 + W2 + N2);

E2 = wvt.Apm_TE_factor.*( abs(wvt.Ap).^2 + abs(wvt.Am).^2 )  + wvt.A0_TE_factor.*( abs(wvt.A0).^2);

sum(E(:))-sum(E2(:))