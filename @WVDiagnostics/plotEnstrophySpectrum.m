function fig = plotEnstrophySpectrum(self,options)
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
% TZ_A0_j_kl = wvt.A0_TZ_factor .* abs(wvt.A0).^2;
prefactor = wvt.h_0/2; prefactor(1) = wvt.Lz/2;
TZ_A0_j_kl = prefactor.*abs(wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.qgpv))).^2;
TZ_APV_j_kl = prefactor.*abs(wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.apv))).^2;
TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);
TZ_APV_j_kR = wvt.transformToRadialWavenumber(TZ_APV_j_kl);
TZ_A0_kR = sum(TZ_A0_j_kR,1);
TZ_A0_j = sum(TZ_A0_j_kR,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% enstrophy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create radial wavelength vector
radialWavelength = 2*pi./wvt.kRadial/1000;
radialWavelength(1) = 2*radialWavelength(2);

% create figure
fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
tl = tiledlayout(2,2,TileSpacing="compact");
title(tl,'Enstrophy Spectrum')

% wave enstrophy???
val = log10((TZ_APV_j_kR).');
axIGW = nexttile;
pcolor(radialWavelength,wvt.j,val.'), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
title('apv')
colormap(axIGW, self.cmocean('dense'));
text(radialWavelength(1),max(wvt.j)*1.05,'MDA','FontWeight','bold')
line([radialWavelength(2),radialWavelength(2)],[min(wvt.j),max(wvt.j)],'Color','k','LineWidth',1.5)
clim([-17 -8])

% plot the geostrophic enstrophy
val = log10((TZ_A0_j_kR).');
axGEO = nexttile;
pcolor(radialWavelength,wvt.j,val.'), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
title('qgpv')
colormap(axGEO, self.cmocean('dense'));
text(radialWavelength(1),max(wvt.j)*1.05,'MDA','FontWeight','bold')
line([radialWavelength(2),radialWavelength(2)],[min(wvt.j),max(wvt.j)],'Color','k','LineWidth',1.5)
clim([-17 -8])

self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

% plot vertical mode spectrum
axJ = nexttile;
plot(wvt.j,TZ_A0_j)
yscale('log')
axis tight
ylabel('enstrophy (m s^{-2})');
xlabel('vertical mode j');
title('Vertical Mode Spectrum')

% plot horizontal wavenumber spectrum
axK = nexttile;
plot(radialWavelength,TZ_A0_kR)
set(gca,'XDir','reverse')
xscale('log'); yscale('log')
axis tight
title('Radial Wavenumber Spectrum')
% ylabel('enstrophy (m s^{-2})');
xlabel('wavelength (km)')

% match limits
xlimK = get(axK,'xlim');
ylimK = get(axK,'ylim');
ylimJ = get(axJ,'ylim');
% set(axIGW,'xlim',xlimK);
set(axGEO,'xlim',xlimK);
set(axJ,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])
set(axK,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])

% colorbar
cb = colorbar(axGEO);
cb.Layout.Tile = 'south';
cb.Label.String = "log10(m s^{-2})";

end
