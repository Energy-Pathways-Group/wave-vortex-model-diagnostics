function [inertial_fluxes,t] = quadraticEnergyTriadFluxesOverTime(self,options)
% Compute inertial fluxes over time
%
% Returns the energy fluxes due to inertial interactions for each reservoir as a function of time.
%
% - Topic: Fluxes over time, [t 1]
% - Declaration: inertial_fluxes = quadraticEnergyTriadFluxesOverTime(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter options.triadComponents: vector of TriadFlowComponent objects (default: [geostrophic_mda, wave])
% - Returns inertial_fluxes: struct array with fluxes over time
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
inertial_fluxes = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end