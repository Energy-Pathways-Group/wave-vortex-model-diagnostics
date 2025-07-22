function S_f = spectrumWithFgTransform(self,f)
arguments
    self WVDiagnostics
    f
end
wvt = self.wvt;
prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor.*abs(f_bar).^2;
end