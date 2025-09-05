function fig = plotEnergyFluxOverTime(self,options)
% Plot forcing flux for each reservoir over time
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
    options.approximation {mustBeMember(options.approximation,{'quadratic','exact'})} = 'exact'
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v) v;
end
if options.approximation == "exact"
    [forcing_fluxes, t] = self.exactEnergyFluxesOverTime(timeIndices=options.timeIndices);

    fig = figure(Visible=options.visible);
    tl = tiledlayout(1,1,TileSpacing="compact");
    for iForce = 1:length(forcing_fluxes)
        plot(t/self.tscale,options.filter(forcing_fluxes(iForce).te/self.flux_scale)), hold on
    end
    legend(forcing_fluxes.fancyName)

    xlabel("time (" + self.tscale_units + ")")
    ylabel("flux (" + self.flux_scale_units + ")")
    xlim([min(t) max(t)]/self.tscale);
else
    [forcing_fluxes, t] = self.quadraticEnergyFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

    fig = figure(Visible=options.visible);
    tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
    for iReservoir = 1:length(options.energyReservoirs)
        nexttile(tl);
        for iForce = 1:length(forcing_fluxes)
            plot(t/self.tscale,options.filter(forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/self.flux_scale)), hold on
        end
        legend(forcing_fluxes.fancyName)

        fancyName = options.energyReservoirs(iReservoir).fancyName;
        xlabel("time (" + self.tscale_units + ")")
        ylabel("flux into " + fancyName + " (" + self.flux_scale_units + ")")
        xlim([min(t) max(t)]/self.tscale);
    end
end
end