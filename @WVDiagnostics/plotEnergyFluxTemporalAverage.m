function fig = plotEnergyFluxTemporalAverage(self,options)
% Plot the wave/geostrophic energy spectra at a given time
%
% Makes a nice multiplanel plot of the wave and geostrophic spectra at a
% given time.
%
% - Topic: Figures (over time)
% - Declaration: fig = plotEnergySpectrum(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter options.triadComponents: vector of TriadFlowComponent objects (default: [geostrophic_mda, wave])
% - Parameter options.timeIndices: indices for time averaging (default: Inf)
% - Parameter options.visible: figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.energyReservoir = EnergyReservoir.total;
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
    options.axes {mustBeMember(options.axes,{'jk','j','k'})} = 'jk'
    options.shouldOverlayWaveFrequencies = true
    options.colormap = WVDiagnostics.crameri('bam')
    options.visible = "on"
    options.overSaturationFactor = 10;
end

wvt = self.wvt;
forcing_fluxes = self.forcingFluxesTemporalAverage(energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
inertial_fluxes = self.inertialFluxesTemporalAverage(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
fluxes = cat(2,forcing_fluxes,inertial_fluxes);

if options.axes == "jk"
    colorLimits = max(arrayfun( @(v) max(abs(v.(options.energyReservoir.name)(:))), fluxes))*[-1 1]/options.overSaturationFactor;
end

fig = figure;
tl = tiledlayout(fig,"flow",TileSpacing='tight');

for iComponent = 1:length(fluxes)
    val = fluxes(iComponent).(options.energyReservoir.name);
    ax = nexttile(tl);
    switch options.axes
        case "jk"
            pcolor(ax,wvt.kRadial,wvt.j,val), shading flat
            colormap(ax, options.colormap)
            self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
            if options.shouldOverlayWaveFrequencies
                self.overlayFrequencyContours;
            end
            clim(ax,colorLimits);
            
        case "j"
            plot(wvt.j,zeros(size(wvt.j)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,wvt.j,sum(val,2))
        case "k"
            plot(wvt.kRadial,zeros(size(wvt.kRadial)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,wvt.kRadial,sum(val,1))
            self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
    end
    title(ax,fluxes(iComponent).fancyName)
end

switch options.axes
    case "jk"
        % jkR plot labels
        xlabel(tl,'wavelength (km)')
        ylabel(tl,"vertical mode")
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = 'energy flux (m^3 s^{-3})';
    case "j"
        % j plot options
        ylabel(tl,'energy flux (m^3 s^{-3})')
        xlabel(tl, 'vertical mode')
    case "k"
        % kR plot options
        ylabel(tl,'energy flux (m^3 s^{-3})')
        xlabel(tl,'wavelength (km)')
end

% remove redundant labels
switch options.axes
    case "jk"
        for ii=1:tilenum(ax) %prod(tl_jkR.GridSize)
            if ii <= prod(tl.GridSize)-tl.GridSize(2)
                nexttile(ii)
                xticklabels([])
            end
            if ~(mod(ii,tl.GridSize(2))==1)
                nexttile(ii)
                yticklabels([])
            end
        end
end
if isinf(options.timeIndices)
    options.timeIndices = 1:length(self.t_diag);
end
minDay = string(round(self.t_diag(min(options.timeIndices))/86400));
maxDay = string(round(self.t_diag(max(options.timeIndices))/86400));
title(tl,"energy flux into " + options.energyReservoir.fancyName + " energy, day " + minDay + "-" + maxDay)

end