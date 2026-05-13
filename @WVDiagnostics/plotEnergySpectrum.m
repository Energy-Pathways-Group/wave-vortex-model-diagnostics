function [fig, spectralSlopes] = plotEnergySpectrum(self,options)
% Plot the wave/geostrophic energy spectra at a given time.
%
% Plot the wave/geostrophic energy spectra at a given time
% Makes a nice multiplanel plot of the wave and geostrophic spectra at a
% given time.
%
% - Topic: Figures — Model Snapshot
% - Declaration: [fig, spectralSlopes] = plotEnergySpectrum(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter iTime: time index in model output file
% - Parameter visible: (optional) figure visibility (default: "on")
% - Returns fig: handle to the generated figure
% - Returns spectralSlopes: struct containing spectral slope fits
%   (`IOIGW_kR_slope`, `A0_kR_slope`, `IOIGW_j_slope`, `A0_j_slope`,
%   `IOIGW_jWavenumber_slope`, `A0_jWavenumber_slope`)
%   Fit lines are shown on the 1D spectra when this output is requested.
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

% create radial wavelength vector
radialWavelength = 2*pi./wvt.kRadial/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

% create j vector for log y-axis.
jForLogAxis = wvt.j;
jForLogAxis(1) = 0.75;

% create figure
fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
tl = tiledlayout(2,2,TileSpacing='tight');
title(tl,'Energy Spectrum')

% plot the wave energy
val = log10(TE_Apm_j_kR);
axIGW = nexttile;
pcolor(radialWavelength,jForLogAxis,val), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
set(gca,'YScale','log')
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
colormap(axIGW, self.cmocean('dense'));
ylabel('vertical mode')
title('Internal Gravity Wave')
xlabel('wavelength (km)')
text(radialWavelength(1),max(jForLogAxis),'IO','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)

self.overlayFrequencyContours(frequencies = [1.01 1.05 1.2 1.5 2 4 8 16],textColor = [.5,.5,.5],labelSpacing = 400,lineWidth = 1)

% plot the geostrophic energy
val = log10(TE_A0_j_kR);
axGEO = nexttile;
pcolor(radialWavelength,jForLogAxis,val), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
set(gca,'YScale','log')
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
colormap(axGEO, self.cmocean('dense'));
set(gca,'YTickLabel',[]);
title('Geostrophic')
xlabel('wavelength (km)')
text(radialWavelength(1),max(jForLogAxis),'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)

self.overlayGeostrophicKineticPotentialFractionContours
self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

% plot vertical mode spectrum
axJ = nexttile;
hJ = plot(wvt.j,TE_inertial_j+TE_wave_j+TE_A0_j,wvt.j,TE_A0_j,wvt.j,TE_inertial_j+TE_wave_j, 'LineWidth',2);
hold on
plot(wvt.j,TE_wave_j,'--','Color',linesTemp(3,:), 'LineWidth',2)
plot(wvt.j,TE_inertial_j,':','Color',linesTemp(3,:), 'LineWidth',2)
xscale('log'),yscale('log')
ylabel('energy (m^3 s^{-2})');
xlabel('vertical mode');
axis tight
title('Vertical Mode Spectrum')
legend('Total','Geostrophic','IO+IGW','IGW','IO','Location','southwest')

% plot horizontal wavenumber spectrum
axK = nexttile;
hK = plot(radialWavelength,TE_inertial_kR+TE_wave_kR+TE_A0_kR,radialWavelength,TE_A0_kR,radialWavelength,TE_inertial_kR+TE_wave_kR, 'LineWidth',2);
hold on
set(gca,'XDir','reverse')
xscale('log'); yscale('log')
axis tight
title('Radial Wavenumber Spectrum')
legend('Total','Geostrophic','IO+IGW','Location','southwest')
xlabel('wavelength (km)')
% yticklabels([])

% match limits
xlimK = get(axK,'xlim');
ylimK = get(axK,'ylim');
ylimJ = get(axJ,'ylim');
set(axIGW,'xlim',xlimK);
set(axGEO,'xlim',xlimK);
set(axJ,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])
set(axJ,'xlim',[1 max(self.j)]) % remove if using linear xlim for axJ
set(axK,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])

% colorbar
cb = colorbar(axIGW);
cb.Layout.Tile = 'south';
cb.Label.String = "log10(m^3 s^{-2})";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% wave spectrum power law fits
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spectralSlopes = struct;

% select wavenumber range
k_no_damp = wvt.forcingWithName('adaptive damping').k_no_damp;
kInd = all([wvt.kRadial>(2*pi/wvt.Lx) , wvt.kRadial<k_no_damp],2);
% select vertical mode range
j_no_damp = wvt.forcingWithName('adaptive damping').j_no_damp;
jInd = all([wvt.j>2 , wvt.j<j_no_damp],2);

% radial wavenumber slope
p_IOIGW_kR = polyfit( log(wvt.kRadial(kInd)), log(TE_inertial_kR(kInd)+TE_wave_kR(kInd)), 1 );
spectralSlopes.IOIGW_kR_slope = p_IOIGW_kR(1);

p_A0_kR = polyfit( log(wvt.kRadial(kInd)), log(TE_A0_kR(kInd)), 1 );
spectralSlopes.A0_kR_slope = p_A0_kR(1);

% j mode slope
p_IOIGW_j = polyfit( log(wvt.j(jInd)), log(TE_inertial_j(jInd)+TE_wave_j(jInd)), 1 );
spectralSlopes.IOIGW_j_slope = p_IOIGW_j(1);

p_A0_j = polyfit( log(wvt.j(jInd)), log(TE_A0_j(jInd)), 1 );
spectralSlopes.A0_j_slope = p_A0_j(1);

% j wavenumber slope
p = polyfit( log(self.jWavenumber(jInd)), log(TE_inertial_j(jInd)+TE_wave_j(jInd)), 1 );
spectralSlopes.IOIGW_jWavenumber_slope = p(1);

p = polyfit( log(self.jWavenumber(jInd)), log(TE_A0_j(jInd)), 1 );
spectralSlopes.A0_jWavenumber_slope = p(1);

% add fit lines to 1d spectrum plots
if nargout > 1
    fitLineOffset = 2;
    fitTextOffset = 2;
    fit_IOIGW_kR = fitLineOffset*exp(polyval(p_IOIGW_kR,log(wvt.kRadial(kInd))));
    fit_A0_kR = fitLineOffset*exp(polyval(p_A0_kR,log(wvt.kRadial(kInd))));
    fit_IOIGW_j = fitLineOffset*exp(polyval(p_IOIGW_j,log(wvt.j(jInd))));
    fit_A0_j = fitLineOffset*exp(polyval(p_A0_j,log(wvt.j(jInd))));

    plot(axK,radialWavelength(kInd),fit_IOIGW_kR, ...
        Color=hK(3).Color,LineStyle=hK(3).LineStyle,LineWidth=1,HandleVisibility='off');
    text(axK,radialWavelength(find(kInd,1,'last')),fitTextOffset*fit_IOIGW_kR(end), ...
        sprintf('%.1f',spectralSlopes.IOIGW_kR_slope),Color=hK(3).Color, ...
        HorizontalAlignment='right',VerticalAlignment='bottom',HandleVisibility='off');

    plot(axK,radialWavelength(kInd),fit_A0_kR, ...
        Color=hK(2).Color,LineStyle=hK(2).LineStyle,LineWidth=1,HandleVisibility='off');
    text(axK,radialWavelength(find(kInd,1,'last')),fitTextOffset*fit_A0_kR(end), ...
        sprintf('%.1f',spectralSlopes.A0_kR_slope),Color=hK(2).Color, ...
        HorizontalAlignment='right',VerticalAlignment='bottom',HandleVisibility='off');

    plot(axJ,wvt.j(jInd),fit_IOIGW_j, ...
        Color=hJ(3).Color,LineStyle=hJ(3).LineStyle,LineWidth=1,HandleVisibility='off');
    text(axJ,wvt.j(find(jInd,1,'last')),fitTextOffset*fit_IOIGW_j(end), ...
        sprintf('%.1f',spectralSlopes.IOIGW_j_slope),Color=hJ(3).Color, ...
        HorizontalAlignment='right',VerticalAlignment='bottom',HandleVisibility='off');

    plot(axJ,wvt.j(jInd),fit_A0_j, ...
        Color=hJ(2).Color,LineStyle=hJ(2).LineStyle,LineWidth=1,HandleVisibility='off');
    text(axJ,wvt.j(find(jInd,1,'last')),fitTextOffset*fit_A0_j(end), ...
        sprintf('%.1f',spectralSlopes.A0_j_slope),Color=hJ(2).Color, ...
        HorizontalAlignment='right',VerticalAlignment='bottom',HandleVisibility='off');
end

end
