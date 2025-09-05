basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

% runNumber=1;
% wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");
wvd = WVDiagnostics(basedir + replace(getRunParameters(17),"256","256") + ".nc");

wvt_lowres = wvd.wvt;
wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();

% wvt = WVTransform.waveVortexTransformFromFile("run18_icR_iner07_tide014_lat32_geo0065_N0052_boussinesq_res512-wvt-aa.nc");
dimensionNames = ["j","kRadial"];
[omegaN,n,hke_jk,pe_jk] = wvt.transformToRadialWavenumber(abs(wvt.Omega),ones(size(wvt.Omega)),wvt.A0_KE_factor,wvt.A0_PE_factor);
omegaJK = (omegaN./n);
wvd.diagfile.addVariable("omega_jk",dimensionNames,omegaJK,isComplex=false);
wvd.diagfile.addVariable("geo_hke_jk",dimensionNames,hke_jk,isComplex=false);
wvd.diagfile.addVariable("geo_pe_jk",dimensionNames,pe_jk,isComplex=false);