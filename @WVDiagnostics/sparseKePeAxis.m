function [kePeAxis,bins_kepe] = sparseKePeAxis(self)
% Sparse Ke Pe Axis.
%
% sparseKePeAxis is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Utilities — Sparse matrices — Axis binning — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
% - Declaration: [kePeAxis,bins_kepe] = sparseKePeAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns kePeAxis: output value `kePeAxis`
% - Returns bins_kepe: output value `bins_kepe`
arguments
    self
end

wvt = self.wvt;

kePeAxis = reshape(self.kePeAxis,1,[]);
kePeFraction = wvt.A0_KE_factor./(wvt.A0_KE_factor+wvt.A0_PE_factor);

mid    = 0.5*(kePeAxis(1:end-1) + kePeAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_kepe = discretize(kePeFraction, edges);
end