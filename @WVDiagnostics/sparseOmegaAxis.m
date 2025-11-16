function [omegaAxis,bins_omega] = sparseOmegaAxis(self)
% Sparse Omega Axis.
%
% sparseOmegaAxis is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Utilities — Sparse matrices — Axis binning — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
% - Declaration: [omegaAxis,bins_omega] = sparseOmegaAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns omegaAxis: output value `omegaAxis`
% - Returns bins_omega: output value `bins_omega`
arguments
    self
end

wvt = self.wvt;

omega = self.omegaAxis;
dlog = (log10(2)-log10(1))/4;
omegaAxis = omega(1)*(10.^(log10(1):dlog:log10(omega(end)/omega(1))));

mid    = 0.5*(omegaAxis(1:end-1) + omegaAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_omega = discretize(wvt.Omega, edges);
end