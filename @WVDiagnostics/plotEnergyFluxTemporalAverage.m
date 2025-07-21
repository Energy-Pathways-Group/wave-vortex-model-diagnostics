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
    options.showForcingFluxes = true;
    options.timeIndices = Inf;
    options.axes {mustBeMember(options.axes,{'jk','j','k'})} = 'jk'
    options.shouldOverlayWaveFrequencies = false %true
    options.shouldOverlayGeostrophicKineticPotentialRatioContours = true
    options.colormap = WVDiagnostics.crameri('-bam')
    options.visible = "on"
    options.overSaturationFactor = 10;
end

n_size = 17;

wvt = self.wvt;
forcing_fluxes = self.forcingFluxesTemporalAverage(energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
inertial_fluxes = self.inertialFluxesTemporalAverage(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
if options.showForcingFluxes
    fluxes = cat(2,forcing_fluxes,inertial_fluxes);
else
    fluxes = inertial_fluxes;
end

if options.axes == "jk"
    colorLimits = max(arrayfun( @(v) max(abs(v.(options.energyReservoir.name)(:))), fluxes))*[-1 1]/options.overSaturationFactor;
    colorLimits = colorLimits/self.flux_scale;
end

% create radial wavelength vector
radialWavelength = 2*pi./wvt.kRadial/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

fig = figure(Visible=options.visible);
tl = tiledlayout(fig,"flow",TileSpacing='tight');

for iComponent = 1:length(fluxes)
    val = fluxes(iComponent).(options.energyReservoir.name)/self.flux_scale;
    ax = nexttile(tl);
    switch options.axes
        case "jk"
            % % % pcolor(ax,options.energyReservoir.kFromKRadial(wvt.kRadial),wvt.j,val), shading flat
            % % % self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
            pcolor(ax,radialWavelength,wvt.j,val), shading flat
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            colormap(ax, options.colormap)
            if options.energyReservoir==EnergyReservoir.total
                text(radialWavelength(1)*.95,0.5,{'MDA','Inertial'},'FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
            elseif options.energyReservoir==EnergyReservoir.geostrophic_mda
                text(radialWavelength(1)*.95,0.5,'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
            elseif options.energyReservoir==EnergyReservoir.wave
                text(radialWavelength(1)*.95,0.5,'Inertial','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
            end
            line([radialWavelength(2),radialWavelength(2)],[min(wvt.j),max(wvt.j)],'Color','k','LineWidth',1)           
            if options.shouldOverlayWaveFrequencies
                self.overlayFrequencyContours;
            end
            if options.shouldOverlayGeostrophicKineticPotentialRatioContours
                % self.overlayGeostrophicKineticPotentialRatioContours;
                self.overlayGeostrophicKineticPotentialFractionContours;
            end
            clim(ax,colorLimits);
            
        case "j"
            plot(wvt.j,zeros(size(wvt.j)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,wvt.j,sum(val,2))
        case "k"
            % % % plot(options.energyReservoir.kFromKRadial(wvt.kRadial),zeros(size(wvt.kRadial)),LineWidth=2,Color=0*[1 1 1]), hold on
            % % % plot(ax,wvt.kRadial,sum(val,1))
            % % % self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
            plot(radialWavelength,zeros(size(radialWavelength)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,radialWavelength,sum(val,1))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')       
    end
    title(ax,fluxes(iComponent).fancyName + " (" + string(sum(val(:))) + " " + self.flux_scale_units + ")" )
end

switch options.axes
    case "jk"
        % jkR plot labels
        xlabel(tl,'wavelength (km)')
        ylabel(tl,"vertical mode")
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = "energy flux (" + self.flux_scale_units + ")";
    case "j"
        % j plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl, 'vertical mode')
    case "k"
        % kR plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl,'wavelength (km)')
end

% remove redundant labels
% switch options.axes
%     case "jk"
%         for ii=1:tilenum(ax) %prod(tl_jkR.GridSize)
%             if ii <= prod(tl.GridSize)-tl.GridSize(2)
%                 nexttile(ii)
%                 xticklabels([])
%             end
%             if ~(mod(ii,tl.GridSize(2))==1)
%                 nexttile(ii)
%                 yticklabels([])
%             end
%         end
% end


switch options.axes
    case "jk"
        for ii=1:tilenum(ax) %prod(tl_jkR.GridSize)
            if ii <= prod(tl.GridSize)-tl.GridSize(2)
                nexttile(ii)
                %xticklabels([])
                % set(gca,'XTickLabels',[],'FontSize',n_size)
                set(gca,'XTickLabels',[])
            end
            if ~(mod(ii,tl.GridSize(2))==1)
                nexttile(ii)
                %yticklabels([])
                % set(gca,'YTickLabels',[],'FontSize',n_size)
                set(gca,'YTickLabels',[])
            end
           if (mod(ii,tl.GridSize(2))==1)
                nexttile(ii)
                set(gca,'YTick',[0 5 10 15])
                % set(gca,'YTickLabels',get(gca,'ytick'),'FontSize',n_size)
                set(gca,'YTickLabels',get(gca,'ytick'))
                % set(gca,'fontname','times')
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