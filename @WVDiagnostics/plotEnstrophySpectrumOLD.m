function fig = plotEnstrophySpectrumOLD(self,options)
% Plot the geostrophic enstrophy spectrum at a given time
%
% Makes a geostrophic enstrophy spectrum at a given time
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

% total enstrophy, A0
TZ_A0_j_kl = wvt.A0_TZ_factor .* abs(wvt.A0).^2;
TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);
TZ_A0_kR = sum(TZ_A0_j_kR,1);
TZ_A0_j = sum(TZ_A0_j_kR,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% enstrophy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create figure
fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
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

self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)

self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

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
self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
xlabel('wavelength (km)')

end
