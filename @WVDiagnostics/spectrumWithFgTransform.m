function S_f = spectrumWithFgTransform(self,f,options)
arguments
    self WVDiagnostics
    f
    options.useExplicitAntialiasedWVT = false
end
if options.useExplicitAntialiasedWVT
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end
prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor.*abs(f_bar).^2;
end