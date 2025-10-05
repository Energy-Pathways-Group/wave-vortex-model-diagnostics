function [j,k,bins_0,bins_pm] = sparseJKAxis(self)
wvt = self.wvt;

k = self.sparseKRadialAxis;
j = self.sparseJWavenumberAxis;

k_edges  = [-Inf; 0.5*(k(1:end-1) + k(2:end)); +Inf];
j_edges  = [-Inf; 0.5*(j(1:end-1) + j(2:end)); +Inf];

[J_ll,K_ll] = ndgrid(j_edges(1:end-1),k_edges(1:end-1));
[J_ur,K_ur] = ndgrid(j_edges(2:end), k_edges(2:end));

jWavenumber = 1./sqrt(wvt.Lr2);
jWavenumber(1) = 0;
J_0 = repmat(jWavenumber,[1 wvt.Nkl]);
Kh = wvt.Kh;
J_pm = 1./sqrt(wvt.g*wvt.h_pm/wvt.f/wvt.f);
J_pm(1,:) = 0; % barotropic mode is a mean?

bins_0 = zeros(wvt.spectralMatrixSize);
bins_pm = zeros(wvt.spectralMatrixSize);
for i=1:length(K_ll(:))
    bins_0(Kh >= K_ll(i) & Kh < K_ur(i) & J_0 >= J_ll(i) & J_0 < J_ur(i)) = i;
    bins_pm(Kh >= K_ll(i) & Kh < K_ur(i) & J_pm >= J_ll(i) & J_pm < J_ur(i)) = i;
end
end