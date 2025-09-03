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
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);
[u,v,ape,apv] = wvt.variableWithName('u','v','ape','apv');
ke = (u.^2 + v.^2)/2;
int_vol(ke + ape) - sum(E(:))
int_vol(ape) - sum(0.5*N2(:))

%%
if isa(wvt.rhoFunction,'chebfun')
    rho_nm = wvt.rhoFunction/wvt.rho0 - 1;
else
    rho_nm = chebfun(wvt.rhoFunction,[min(wvt.z) max(wvt.z)],'splitting','on')/wvt.rho0 - 1;
end
p_nm = - wvt.g * cumsum(rho_nm);
p_nm = p_nm - p_nm(0);
eta_true = wvt.eta_true;
Z = wvt.Z;
int_vol(wvt.g*eta_true.*rho_nm(Z - eta_true) + p_nm(Z) - p_nm(Z - eta_true))

first_term = wvd.crossSpectrumWithGgTransform(wvt.g*eta_true./wvt.N2Function(Z),rho_nm(Z - eta_true));
second_term = wvd.crossSpectrumWithFgTransform(p_nm(Z) - p_nm(Z - eta_true),ones(size(Z)));
spectrum = first_term + second_term;
sum(spectrum(:))

%%
alt_first_term = 0.5*wvd.crossSpectrumWithGgTransform(eta_true.*wvt.N2Function(Z-eta_true)./wvt.N2Function(Z),eta_true);
sum(alt_first_term(:))
%%

prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(p_nm(Z) - p_nm(Z - eta_true)));
S_f = prefactor.*abs(f_bar).^2;
sum(S_f(:))


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