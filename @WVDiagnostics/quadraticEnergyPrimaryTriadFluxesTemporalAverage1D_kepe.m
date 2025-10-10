function [inertial_fluxes_g, kePeAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

% [M_wwg, kePeAxis] = self.quadraticEnergyMirrorTriadFluxes1D_omega(timeIndices=options.timeIndices);
% flux_interp = @(v) diff(cat(1,zeros(1,size(v,2)),interp1(self.omegaAxis,cumsum(v),kePeAxis)));

kePeAxis = self.kePeAxis;

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_g = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.geostrophic_mda,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_g)
    fluxes_g(idx).flux = self.transformToKePeAxis(fluxes_g(idx).te_gmda);
    % fluxes_g(idx).flux = flux_interp(val);
end

inertial_fluxes_g(1).flux = fluxes_g([fluxes_g.name] == "gmda_gmda").flux;
inertial_fluxes_g(1).name = "ggg";
inertial_fluxes_g(1).fancyName = "$[g{\nabla}g]_g$";

inertial_fluxes_g(2).flux = 0*kePeAxis; %fluxes_g([fluxes_g.name] == "gmda_wave").flux + fluxes_g([fluxes_g.name] == "wave_gmda").flux; % + M_ggw;
inertial_fluxes_g(2).name = "ggw";
inertial_fluxes_g(2).fancyName = "$[g{\nabla}w]_g+[w{\nabla}g]_g+\mathcal{M} [g{\nabla}g]_w$ (incomplete!!!)";

inertial_fluxes_g(3).flux = fluxes_g([fluxes_g.name] == "wave_wave").flux;
inertial_fluxes_g(3).name = "tx-wwg";
inertial_fluxes_g(3).fancyName = "$[w{\nabla}w]_g$";

inertial_fluxes_g(4).flux = 0*kePeAxis; %-M_ggw;
inertial_fluxes_g(4).name = "tx-ggw";
inertial_fluxes_g(4).fancyName = "$-\mathcal{M} [g{\nabla}g]_w$ (incomplete!!!)";

end