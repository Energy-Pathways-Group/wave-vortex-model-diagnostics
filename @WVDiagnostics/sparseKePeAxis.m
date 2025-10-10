function [kePeAxis,bins_kepe] = sparseKePeAxis(self)
wvt = self.wvt;

kePeAxis = reshape(self.kePeAxis,1,[]);
kePeFraction = wvt.A0_KE_factor./(wvt.A0_KE_factor+wvt.A0_PE_factor);

mid    = 0.5*(kePeAxis(1:end-1) + kePeAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_kepe = discretize(kePeFraction, edges);
end