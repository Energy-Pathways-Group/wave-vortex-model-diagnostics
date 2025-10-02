function [omegaAxis,bins_omega] = sparseOmegaAxis(self)
wvt = self.wvt;

omega = self.omegaAxis;
dlog = (log10(2)-log10(1))/4;
omegaAxis = omega(1)*(10.^(log10(1):dlog:log10(omega(end)/omega(1))));

mid    = 0.5*(omegaAxis(1:end-1) + omegaAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_omega = discretize(wvt.Omega, edges);
end