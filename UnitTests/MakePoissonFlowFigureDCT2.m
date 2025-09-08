basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=1; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
energy_fluxes = wvd.exactEnergyFluxesTemporalAverage(timeIndices=51:251);


%%
flux = (energy_fluxes(1).te.')/wvd.flux_scale;

% For the DCT2/DST2 we use a half-shift grid
x = wvd.kRadial + (wvd.kRadial(2)-wvd.kRadial(1))/2;
y = wvd.jWavenumber + (wvd.jWavenumber(2)-wvd.jWavenumber(1))/2;

N = length(x);
M = length(y);
dk = 1/(2*N*(x(2)-x(1)));
dl = 1/(2*M*(y(2)-y(1)));

k=dk*(0:(N-1)).';
l=dl*(0:(M-1)).';

DCTx = WVDiagnostics.DCT2(N);
DCTy = WVDiagnostics.DCT2(M);

flux_ky = DCTx*flux;
flux_kl = shiftdim(DCTy*shiftdim(flux_ky,1),1);

[K,L] = ndgrid(k,l);

D = -((2*pi*K).^2 + (2*pi*L).^2);
D(1,1) = Inf;

UFactor = 2*pi*K./D;
VFactor = 2*pi*L./D;

iDCTx = WVDiagnostics.iDCT2(N);
iDSTy = WVDiagnostics.iDST2(M);
iDSTy = circshift(iDSTy,1,2);
iDSTy(:,1) = 0;

V_xl = iDCTx*(VFactor.*flux_kl);
V = shiftdim(iDSTy*shiftdim(V_xl,1),1);

iDCTy = WVDiagnostics.iDCT2(M);
iDSTx = WVDiagnostics.iDST2(N);
iDSTx = circshift(iDSTx,1,2);
iDSTx(:,1) = 0;

U_xl = iDSTx*(UFactor.*flux_kl);
U = shiftdim(iDCTy*shiftdim(U_xl,1),1);

[X,Y] = ndgrid(x,y);
figure
quiver(X,Y,U,V,Color=0*[1 1 1])