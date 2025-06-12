function [fig_energy, fig_enstrophy] = plotEnergySpectrum(self,options)
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
TE_Apm_kR = sum(TE_Apm_j_kR,1);
TE_Apm_j = sum(TE_Apm_j_kR,2);

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


% total enstrophy, A0
TZ_A0_j_kl = wvt.A0_TZ_factor .* abs(wvt.A0).^2;
TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);
TZ_A0_kR = sum(TZ_A0_j_kR,1);
TZ_A0_j = sum(TZ_A0_j_kR,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Combined energy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default line color order, for repeating
linesTemp = lines;

% color for some annotations
textColor = '[.5,.5,.5]';

% create some nice tick labels to show wavelength
num_ticks = 6;
roundToNearest = 5;
ticks_x = logspace(log10(wvt.kRadial(2)),log10(wvt.kRadial(end)),num_ticks);
ticks_x = round(2*pi./(1e3.*ticks_x)/roundToNearest)*roundToNearest;
labels_x = cell(length(ticks_x),1);
for i=1:length(ticks_x)
    labels_x{i} = sprintf('%.0f',ticks_x(i));
end
ticks_x = 2*pi./(1e3*ticks_x);

% create the lines of constant frequency
[omegaN,n] = wvt.transformToRadialWavenumber(abs(wvt.Omega),ones(size(wvt.Omega)));
omegaJK = (omegaN./n)/wvt.f;

% create figure
fig_energy = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
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
xscale('log')
xticks(ticks_x)
xticklabels(labels_x)
set(gca,'YTickLabel',[]);
title('wave')
xlabel('wavelength (km)')
% add frequency contours
set(gca,'layer','top'),
hold on
v = [1.01 1.05 1.2 1.5 2 4 8 16];
[C,h] = contour(wvt.kRadial(2:end),wvt.j',(omegaJK(:,2:end)),v,'LineWidth',1,'Color',textColor);
clabel(C,h,v,'Color',textColor,'LabelSpacing',400)
% % val = log10((TE_Apm_j_kR(:,2:end)));
% % ax = nexttile([1 inertialPlotWidthRatio]);
% % pcolor(2*pi/1000./wvt.kRadial(2:end),wvt.j,val), shading flat,
% % clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
% % % clim([max(var(:))-6 max(var(:))])
% % colormap(ax, self.cmocean('dense'));
% % % xscale('log')
% % set(gca,'xdir','reverse')
% % % set(gca,'XScale','log')
% % % set(gca,'XLim',[1e1,1e3])
% % % xt = xticks;
% % % xtl = compose('%g', round(10.^xt));
% % % set(gca,'XTickLabel',xtl)
% % % xticks(ticks_x)gca
% % % xticklabels(labels_x)
% % set(gca,'YTickLabel',[]);
% % title('wave')
% % xlabel('wavelength (km)')
% % % add frequency contours
% % set(gca,'layer','top'),
% % hold on
% % v = [1.01 1.05 1.2 1.5 2 4 8 16];
% % [C,h] = contour(wvt.kRadial(2:end),wvt.j',(omegaJK(:,2:end)),v,'LineWidth',1,'Color',textColor);
% % clabel(C,h,v,'Color',textColor,'LabelSpacing',400

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
set(gca,'Layer','top','TickLength',[0.015 0.015])
% create some nice tick labels to show deformation radius
yticksTemp = yticks;
ticks_y = sqrt(wvt.Lr2)./1000;
labels_y = cell(length(yticksTemp),1);
for i=1:length(yticksTemp)
    labels_y{i} = sprintf('%0.1f',ticks_y(yticksTemp(i)+1));
end
text(1.25*max(xlim)*ones(size(yticksTemp)),yticksTemp,labels_y,'Color',textColor,'HorizontalAlignment','left')
text(2*max(xlim),1.1*max(ylim),'L_r (km)','Color',textColor,'HorizontalAlignment','right')

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

% save figure
% exportgraphics(gcf,figFolder + "energy_spectrum.png")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% enstrophy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create figure
fig_enstrophy = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
inertialPlotWidthRatio = 4;
tl = tiledlayout(2,2*inertialPlotWidthRatio+1,TileSpacing="compact");
title(tl,'Enstrophy Spectrum')

% wave enstrophy???

% plot the geostrophic enstrophy
val = log10((TZ_A0_j_kR).');
ax = nexttile(inertialPlotWidthRatio+2, [1 inertialPlotWidthRatio]);
pcolor(wvt.kRadial,wvt.j,val.'), shading flat,
% clim([max(var(:))-6 max(var(:))])
colormap(ax, self.cmocean('dense'));
xscale('log')
xticks(ticks_x)
xticklabels(labels_x)
set(gca,'YTickLabel',[]);
title('geostrophic')
xlabel('wavelength (km)')
set(gca,'layer','top'), hold on
% create some nice tick labels to show deformation radius
yticksTemp = yticks;
ticks_y = sqrt(wvt.Lr2)./1000;
labels_y = cell(length(yticksTemp),1);
for i=1:length(yticksTemp)
    labels_y{i} = sprintf('%0.1f',ticks_y(yticksTemp(i)+1));
end
text(1.25*max(xlim)*ones(size(yticksTemp)),yticksTemp,labels_y,'Color',textColor,'HorizontalAlignment','left')
text(2*max(xlim),1.1*max(ylim),'L_r (km)','Color',textColor,'HorizontalAlignment','right')

% plot vertical mode spectrum
ax = nexttile(2*inertialPlotWidthRatio+3,[1 inertialPlotWidthRatio]);
plot(wvt.j,TZ_A0_j)
yscale('log')
ylabel('enstrophy (m s^{-2})');
xlabel('vertical mode j');
title('Vertical Mode Spectrum')

% plot horizontal wavenumber spectrum
ax = nexttile(3*inertialPlotWidthRatio+3,[1 inertialPlotWidthRatio]);
plot(wvt.kRadial,TZ_A0_kR)
xscale('log'); yscale('log')
title('Radial Wavenumber Spectrum')
ylabel('enstrophy (m s^{-2})');
xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)')

end
