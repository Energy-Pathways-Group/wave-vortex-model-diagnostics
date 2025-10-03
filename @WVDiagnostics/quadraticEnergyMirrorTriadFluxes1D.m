function [M_wwg, M_ggw, kp] = quadraticEnergyMirrorTriadFluxes1D(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
% The mirror fluxes of the the [g{\nabla}g]_w and [w{\nabla}w]_g triad
% components are computed on a custom sparse pseudo-radial wavelength grid

val= self.diagfile.readVariables('pi_w_wwg_kp');
M_wwg= cat(1,zeros(1,size(val,2)),diff(val));

val= self.diagfile.readVariables('pi_g_ggw_kp');
M_ggw= cat(1,zeros(1,size(val,2)),diff(val));

kp =  reshape(self.diagfile.readVariables('kp'),[],1);

if isinf(options.timeIndices)
    M_wwg = mean(M_wwg,2);
    M_ggw = mean(M_ggw,2);
else
    M_wwg = mean(M_wwg(:,options.timeIndices),2);
    M_ggw = mean(M_ggw(:,options.timeIndices),2);
end

end