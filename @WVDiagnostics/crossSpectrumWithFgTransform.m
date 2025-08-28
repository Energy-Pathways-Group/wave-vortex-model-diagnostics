function S_f = crossSpectrumWithFgTransform(self,phi,gamma)
arguments
    self WVDiagnostics
    phi
    gamma
end
wvt = self.wvt;
% dp=WVGeometryDoublyPeriodic([wvt.Lx wvt.Ly],[wvt.Nx wvt.Ny],Nz=wvt.Nz,conjugateDimension=wvt.conjugateDimension,shouldAntialias=false);

prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

phi_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(phi));
gamma_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(gamma));
S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
end