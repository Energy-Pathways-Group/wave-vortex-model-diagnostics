function [M_wwg, omegaAxis] = quadraticEnergyMirrorTriadFluxes1D_omega(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
% The mirror fluxes of the the [w{\nabla}w]_g triad
% components are computed on a custom sparse omega grid

val= self.diagfile.readVariables('pi_w_wwg_omega');
M_wwg= cat(1,zeros(1,size(val,2)),diff(val));

omegaAxis =  reshape(self.diagfile.readVariables('omegaAxis'),[],1);

if isinf(options.timeIndices)
    M_wwg = mean(M_wwg,2);
else
    M_wwg = mean(M_wwg(:,options.timeIndices),2);
end

end