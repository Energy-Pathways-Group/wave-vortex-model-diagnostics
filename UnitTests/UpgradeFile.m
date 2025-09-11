basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=18;
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");
% wvd = WVDiagnostics(basedir + replace(getRunParameters(17),"256","256") + ".nc");

%%
% wvt_lowres = wvd.wvt;
% wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
wvt = WVTransform.waveVortexTransformFromFile(basedir+"run18_icR_iner07_tide014_lat32_geo0065_N0052_boussinesq_res512-wvt-aa.nc");

if isa(wvt,"WVTransformHydrostatic")
    h_pm = repmat(wvt.h_pm,[1 wvt.Nkl]);
else
    h_pm = wvt.h_pm;
end

% wvt = WVTransform.waveVortexTransformFromFile("run18_icR_iner07_tide014_lat32_geo0065_N0052_boussinesq_res512-wvt-aa.nc");
dimensionNames = ["j","kRadial"];
[n,hN] = wvt.transformToRadialWavenumber(ones(size(wvt.Omega)),h_pm);
h_kj = (hN./n);
Lr2_pm = wvt.g * h_kj / (wvt.f*wvt.f);
wvd.diagfile.addVariable("Lr2_pm",dimensionNames,Lr2_pm,isComplex=false);