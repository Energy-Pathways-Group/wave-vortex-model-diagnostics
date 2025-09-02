function enstrophy_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
% Compute spatial-temporal average of the exact enstrophy fluxes
%
% Returns the spatial-temporal average of the exact enstrophy fluxes from external forcing
%
% - Topic: Flux averages, scalar
% - Declaration: forcing_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
% - Parameter options.timeIndices: indices for time averaging (default: Inf)
% - Returns forcing_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

enstrophy_fluxes = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
for iForce = 1:length(enstrophy_fluxes)
    enstrophy_fluxes(iForce).Z0 = mean(enstrophy_fluxes(iForce).Z0);
end
end