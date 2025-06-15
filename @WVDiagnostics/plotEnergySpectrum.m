function fig = plotEnergySpectrum(self,options)
% Plot the wave/geostrophic energy spectra at a given time
%
% Makes a nice multiplanel plot of the wave and geostrophic spectra at a
% given time.
%
% - Topic: Figures (over time)
% - Declaration: fig = plotEnergySpectrum(self,options)
% - Parameter options.iTime: time index in model output file
% - Parameter options.visible: figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.iTime
    options.visible = "on"
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% compute energy and enstrophy spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total energy, A0
TE_A0_j_kl = wvt.A0_TE_factor .* abs(wvt.A0).^2; % m^2/s^3
TE_A0_j_kR = wvt.transformToRadialWavenumber(TE_A0_j_kl);
TE_A0_kR = sum(TE_A0_j_kR,1);
TE_A0_j = sum(TE_A0_j_kR,2);

% total energy, Apm
TE_Apm_j_kl = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2); % m^2/s^3
TE_Apm_j_kR = wvt.transformToRadialWavenumber(TE_Apm_j_kl);

% total energy, inertial
maskApInertial = wvt.inertialComponent.maskAp; % m^2/s^3
maskAmInertial = wvt.inertialComponent.maskAm;
TE_inertial_j_kl = wvt.Apm_TE_factor .* (abs(maskApInertial.*wvt.Ap).^2 + abs(maskAmInertial.*wvt.Am).^2);
TE_inertial_j_kR = wvt.transformToRadialWavenumber(TE_inertial_j_kl);
TE_inertial_kR = sum(TE_inertial_j_kR,1);
TE_inertial_j = sum(TE_inertial_j_kR,2);

% total energy, wave
maskApWave = wvt.waveComponent.maskAp; % m^2/s^3
maskAmWave = wvt.waveComponent.maskAm;
TE_wave_j_kl = wvt.Apm_TE_factor .* (abs(maskApWave.*wvt.Ap).^2 + abs(maskAmWave.*wvt.Am).^2);
TE_wave_j_kR = wvt.transformToRadialWavenumber(TE_wave_j_kl);
TE_wave_kR = sum(TE_wave_j_kR,1);
TE_wave_j = sum(TE_wave_j_kR,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Combined energy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default line color order, for repeating
linesTemp = lines;

% create figure
fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
inertialPlotWidthRatio = 4;
tl = tiledlayout(2,2*inertialPlotWidthRatio+1,TileSpacing="tight");
title(tl,'Energy Spectrum')

% plot the inertial energy
val = log10(repmat(TE_Apm_j_kR(:,1),[1 2]));
ax = nexttile;
pcolor([0;1],wvt.j,val), shading flat,
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
% clim([max(var(:))-6 max(var(:))])
colormap(ax, self.cmocean('dense'));
set(gca,'XTickLabel',[]);
ylabel('mode')
title('inertial')
set(gca,'Layer','top','TickLength',[0.015 0.015])

% plot the wave energy
val = log10((TE_Apm_j_kR(:,2:end)));
ax = nexttile([1 inertialPlotWidthRatio]);
pcolor(wvt.kRadial(2:end),wvt.j,val), shading flat,
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
% clim([max(var(:))-6 max(var(:))])
colormap(ax, self.cmocean('dense'));

self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)

set(gca,'YTickLabel',[]);
title('wave')
xlabel('wavelength (km)')

self.overlayFrequencyContours(frequencies = [1.01 1.05 1.2 1.5 2 4 8 16],textColor = [.5,.5,.5],labelSpacing = 400,lineWidth = 1)

% plot the geostrophic energy
val = log10((TE_A0_j_kR).');
ax = nexttile([1 inertialPlotWidthRatio]);
pcolor(wvt.kRadial,wvt.j,val.'), shading flat,
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
% clim([max(var(:))-6 max(var(:))])
colormap(ax, self.cmocean('dense'));
xscale('log')
xticks(ticks_x)
xticklabels(labels_x)
set(gca,'YTickLabel',[]);
title('geostrophic')
xlabel('wavelength (km)')

self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

% plot vertical mode spectrum
ax = nexttile(2*inertialPlotWidthRatio+3,[1 inertialPlotWidthRatio]);
plot(wvt.j,TE_inertial_j+TE_wave_j+TE_A0_j,wvt.j,TE_A0_j,wvt.j,TE_inertial_j+TE_wave_j)
hold on
plot(wvt.j,TE_wave_j,'--','Color',linesTemp(3,:))
plot(wvt.j,TE_inertial_j,':','Color',linesTemp(3,:))
yscale('log')
ylabel('energy (m^3 s^{-2})');
xlabel('vertical mode j');
title('Vertical Mode Spectrum')
legend('Total','Geostrophic','IO+IGW','IGW','IO','Location','southwest')

% plot horizontal wavenumber spectrum
ax = nexttile([1 inertialPlotWidthRatio]);
plot(wvt.kRadial,TE_inertial_kR+TE_wave_kR+TE_A0_kR,wvt.kRadial,TE_A0_kR,wvt.kRadial,TE_inertial_kR+TE_wave_kR)
xscale('log'); yscale('log')
title('Radial Wavenumber Spectrum')
legend('Total','Geostrophic','IO+IGW','Location','southwest')
xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)')
yticklabels([])
ylabel('energy (m^3 s^{-2})');

end
