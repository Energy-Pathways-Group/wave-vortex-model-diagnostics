function inertial_fluxes = quadraticEnergyTriadFluxesSpatialTemporalAverage(self,options)
% Compute spatial-temporal average of inertial fluxes
%
% Returns the spatial-temporal average of energy fluxes due to inertial interactions for each reservoir.
%
% - Topic: Flux averages, scalar
% - Declaration: inertial_fluxes = quadraticEnergyTriadFluxesSpatialTemporalAverage(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter options.timeIndices: indices for time averaging (default: Inf)
% - Parameter options.triadComponents: vector of TriadFlowComponent objects (default: [geostrophic_mda, wave])
% - Returns inertial_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
end

if isinf(options.timeIndices)
    filter_space = @(v) sum(sum(mean(v,3),1),2);
else
    filter_space = @(v) sum(sum(mean(v(:,:,options.timeIndices),3),1),2);
end
inertial_fluxes = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
end