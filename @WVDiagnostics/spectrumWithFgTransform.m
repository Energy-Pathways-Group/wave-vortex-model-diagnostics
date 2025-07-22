function S_f = spectrumWithFgTransform(self,f)
arguments
    self WVDiagnostics
    f
end
wvt = self.wvt;
prefactor = wvt.h_0; prefactor(1) = wvt.Lz;
f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor.*abs(f_bar).^2;
end