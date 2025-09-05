function fig = plotEnstrophyTriadFluxOverTime(self,options)
% Plot forcing flux for each reservoir over time
%
% filter=@(v,t) movmean(v,51)
%
% Plots the energy flux into each reservoir from external forcing as a function of time.
%
% - Topic: Figures (over time)
% - Declaration: fig = plotForcingFluxOverTime(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter options.visible: figure visibility (default: "on")
% - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v,t) v;
end
[forcing_fluxes, t] = self.quadraticEnstrophyTriadFluxesOverTime(timeIndices=options.timeIndices,triadComponents=options.triadComponents);

fig = figure(Visible=options.visible);
tl = tiledlayout(1,1,TileSpacing="compact");

for iForce = 1:length(forcing_fluxes)
    plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale,t)), hold on
end
legend(forcing_fluxes.fancyName)

xlabel("time (" + self.tscale_units + ")")
ylabel("enstrophy flux (" + self.z_flux_scale_units + ")")
xlim([min(t) max(t)]/self.tscale);

end