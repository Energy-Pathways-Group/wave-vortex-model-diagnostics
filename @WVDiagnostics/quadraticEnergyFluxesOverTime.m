function [forcing_fluxes,t] = quadraticEnergyFluxesOverTime(self,options)
% Compute forcing fluxes over time
%
% Returns the energy fluxes from external forcing for each reservoir as a function of time.
%
% - Topic: Fluxes over time, [t 1]
% - Declaration: forcing_fluxes = quadraticEnergyFluxesOverTime(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter options.timeIndices: indices for time selection (default: Inf)
% - Returns forcing_fluxes: struct array with fluxes over time
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
forcing_fluxes = self.quadraticEnergyFluxes(energyReservoirs=options.energyReservoirs);
exact_forcing_fluxes = self.exactForcingFluxesOverTime();
for iForce=1:length(forcing_fluxes)
    forcing_fluxes(iForce).te = reshape(exact_forcing_fluxes(iForce).te,1,1,[]);
end

forcing_fluxes = self.filterFluxesForReservoir(forcing_fluxes,filter=filter_space);
t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end