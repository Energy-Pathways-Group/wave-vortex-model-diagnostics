function fig = plotFluidStateMultipanel(self,options)
% Plot multipanel summary of fluid state and spectra
%
% Create a compact multipanel figure showing horizontal (x-y) maps and
% vertical (x-z) sections of vertical vorticity for the total flow, the
% wave component, and the geostrophic component at a specified model time.
% Optionally includes log-energy spectra with KE/PE and frequency/wavelength
% contours. Axes are annotated in kilometers and depth in kilometers; color
% limits and annotation ticks are handled internally.
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotFluidStateMultipanel(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter visible: (optional) input argument `visible` (default: "on")
% - Parameter iTime: time-related parameter `iTime`
% - Parameter title: input argument `title`
% - Parameter shouldShowEnergySpectra: (optional) input argument `shouldShowEnergySpectra` (default: true)
% - Parameter shouldShowTotalFields: (optional) input argument `shouldShowTotalFields` (default: false)
% - Parameter figureHandle: input argument `figureHandle`
% - Parameter wavelengths: (optional) input argument `wavelengths` (default: [1,2,5,10,20,50,100,200,500])
% - Parameter wavelengthColor: (optional) input argument `wavelengthColor` (default: [.5,.5,.5])
% - Parameter frequencies: (optional) input argument `frequencies` (default: [1.01 1.05 1.2 2 4 8 16])
% - Parameter frequencyColor: (optional) input argument `frequencyColor` (default: [.7,.7,.7])
% - Parameter keFractions: (optional) input argument `keFractions` (default: [.01,.1,.25,.5,.75,.9,.99])
% - Parameter keFractionColor: (optional) input argument `keFractionColor` (default: [.7,.7,.7])
% - Parameter labelSpacing: (optional) input argument `labelSpacing` (default: 1000)
% - Parameter lineWidth: (optional) input argument `lineWidth` (default: 1)
% - Returns fig: Figure handle for the generated plot
arguments
    self WVDiagnostics
    options.visible = "on"
    options.iTime
    options.title
    options.shouldShowEnergySpectra = true
    options.shouldShowTotalFields = false
    options.figureHandle
    options.wavelengths = [1,2,5,10,20,50,100,200,500];
    options.wavelengthColor = [.5,.5,.5];
    options.frequencies = [1.01 1.05 1.2 2 4 8 16];
    options.frequencyColor = [.7,.7,.7];
    options.keFractions = [.01,.1,.25,.5,.75,.9,.99];
    options.keFractionColor = [.7,.7,.7];
    options.labelSpacing = 1000;
    options.lineWidth = 1;
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

if ~isfield(options,"title")
    options.title = sprintf('%d days',round(wvt.t/86400));
end

cmDivRWB = self.cmocean('balance'); % diverging positive-negative
cmLinMono = self.cmocean('dense') ; % linear monochrome for data>0

% set limits
zeta_limits = [-0.2 0.2];
energy_limits = [-8 0];
enstrophy_limits = [-16 -9];

% We will pretend the "0" wavenumber is actually evenly spaced
% from the nearest two wavenumbers
kPseudoLocation = wvt.kRadial;
kPseudoLocation(1) = exp(-log(kPseudoLocation(3)) + 2*log(kPseudoLocation(2)));
kModePseudoLocation = wvt.kRadial/wvt.dk;
kModePseudoLocation(1) = exp(-log(kModePseudoLocation(3)) + 2*log(kModePseudoLocation(2)));
jwnPseudoLocation = self.jWavenumber(wvt.j+1);
jwnPseudoLocation(1) = exp(-log(jwnPseudoLocation(3)) + 2*log(jwnPseudoLocation(2)));
jPseudoLocation = wvt.j;
jPseudoLocation(1) = exp(-log(jPseudoLocation(3)) + 2*log(jPseudoLocation(2)));
[KPseudoLocation,JWNPseudoLocation] = ndgrid(kPseudoLocation,jwnPseudoLocation);
[KModePseudoLocation,JPseudoLocation] = ndgrid(kModePseudoLocation,jPseudoLocation);
KPseudoRadial = sqrt(JWNPseudoLocation.^2 + KPseudoLocation.^2);

% NOTE: these are from wvd, including the anti-aliased modes.
% For interpolation of hke/omega from WVD to work correctly we need to
% repeat the first entry, but properly back at zero
kPseudoLocationWVD = self.kRadial;
kPseudoLocationWVD(1) = exp(-log(kPseudoLocationWVD(3)) + 2*log(kPseudoLocationWVD(2)));
jPseudoLocationWVD = self.j;
jPseudoLocationWVD(1) = exp(-log(jPseudoLocationWVD(3)) + 2*log(jPseudoLocationWVD(2)));
kPaddedWVD = cat(1,0,kPseudoLocationWVD);
jPaddedWVD = cat(1,0,jPseudoLocationWVD);
[KPaddedWVD,JPaddedWVD] = ndgrid(kPaddedWVD,jPaddedWVD);

% % For interpolation to work correctly we need to repeat the
% % first entry, but properly back at zero
% kPadded = cat(1,0,kPseudoLocation);
% kModePadded = cat(1,0,kModePseudoLocation);
% jwnPadded= cat(1,0,jwnPseudoLocation);
% jPadded= cat(1,0,jPseudoLocation);
% [KModePadded,JPadded] = ndgrid(kModePadded,jPadded); % mode number grids
% [KPadded,JWNPadded] = ndgrid(kPadded,jwnPadded); % wavenumber grids

% color for some annotations
textColor = '[.5,.5,.5]';

% location for x-z section
% iY = round(wvt.Nx/2);
iY = 1;


nColumns = 2;
if options.shouldShowEnergySpectra
    nColumns = nColumns + 1;
end
if nColumns == 3
    figPos = [50 50 900 600];
else
    figPos = [50 50 600 615];
end

if ~isfield(options,"figureHandle")
    fig = figure(Units='points',Position=figPos,Visible = options.visible);
    set(gcf,'PaperPositionMode','auto')
else
    fig = options.figureHandle;
    clf(options.figureHandle)
    set(0, 'currentfigure', options.figureHandle);
end

tl = tiledlayout(1,nColumns,TileSpacing="tight");
title(tl, options.title, 'Interpreter', 'none')

% compute some quantities
TE_A0_j_kl = wvt.A0_TE_factor .* abs(wvt.A0).^2; % m^2/s^3
TE_A0_j_kR = wvt.transformToRadialWavenumber(TE_A0_j_kl);
TE_Apm_j_kl = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2); % m^2/s^3
TE_Apm_j_kR = wvt.transformToRadialWavenumber(TE_Apm_j_kl);
TZ_A0_j_kl = wvt.A0_TZ_factor .* (wvt.A0.*conj(wvt.A0));
TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);

% wave and geosgrophic vorticity v_x - u_y
zeta_z_g = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);  % geostrophic
zeta_z_w = wvt.diffX(wvt.v_w) - wvt.diffY(wvt.u_w);  % wave

% nested tiled layout allows common colorbar for subset of axes.
tl_inner = tiledlayout(tl,2,2,TileSpacing='tight');
tl_inner.Layout.TileSpan = [1,2];

% geostrophic surface vorticity
ax = nexttile(tl_inner,1);
val = zeta_z_g(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:');hold off; % add line for x-z section
title("geostrophic vorticity \zeta_g")
axis square
xticklabels([])
ylabel('y-distance (km)')
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
colormap(ax, cmDivRWB);
clim(ax, zeta_limits);

% geostrophic vorticity section
ax = nexttile(tl_inner,3);
val = squeeze(zeta_z_g(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
xlabel('x-distance (km)')
ylabel('Depth (km)')
axis square
colormap(ax, cmDivRWB);
ylabel('Depth (km)')
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);

% wave surface vorticity
ax = nexttile(tl_inner,2);
val = zeta_z_w(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:');hold off; % add line for x-z section
title("wave vorticity \zeta_w")
axis square
colormap(ax, cmDivRWB);
yticklabels([])
xticklabels([])
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);

% wave vorticity section
ax = nexttile(tl_inner,4);
val = squeeze(zeta_z_w(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
axis square
colormap(ax, cmDivRWB);
xlabel('x-distance (km)')
yticklabels([])
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.String = "$(f)$"; %"$\zeta$/f";
cb.Label.Interpreter = 'latex';

if options.shouldShowEnergySpectra

    % nested tiled layout allows common colorbar for subset of axes.
    tl_inner = tiledlayout(tl,2,1,TileSpacing='tight');
    tl_inner.Layout.Tile = 3;

    % geostrophic energy spectrum
    ax = nexttile(tl_inner,1);
    val = log10(TE_A0_j_kR);
    pcolor(ax,2*pi./kPseudoLocation/1000,jPseudoLocation,val), shading flat,
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xticklabels([])
    title('geostrophic energy spectrum')
    axis square
    clim(ax,energy_limits);
    ylabel("Vertical mode number")
    set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
    colormap(ax, cmLinMono);
    set(gca,'layer','top'),
    hold on
    % add ke:pe ratio contours. flipud/fliplr gives nicer clabel placement.
    % Note: have to pad quantities to work right with our log-log axes. And here
    % have to "double pad" to work right with pcolor's box shift and
    % reversed axis direction. 
    fraction = self.geo_hke_jk./(self.geo_hke_jk+self.geo_pe_jk);
    fractionPadded = cat(1,[fraction(1,:);fraction(1,:)],fraction(1:end-1,:));
    fractionPadded = cat(2,[fractionPadded(:,1) fractionPadded(:,1)],fractionPadded(:,1:end-1));
    fractionJK = interpn(KPaddedWVD,JPaddedWVD,fractionPadded', KPseudoLocation,JPseudoLocation,'linear');
    [C,h] = contour(ax,flipud(2*pi./kPseudoLocation/1000),jPseudoLocation,fliplr(fractionJK'),options.keFractions,'LineWidth',options.lineWidth,'Color',options.keFractionColor, DisplayName="KE/(KE+PE)", HandleVisibility='off');
    clabel(C,h,options.keFractions,'Color',options.keFractionColor,'LabelSpacing',options.labelSpacing)
    % add pseudoWavelength. flipud/fliplr gives nicer clabel placement.
    [C,h] = contour(ax,flipud(2*pi./kPseudoLocation/1000),jPseudoLocation,fliplr(2*pi./KPseudoRadial'/1000),options.wavelengths,'LineWidth',options.lineWidth,'Color',options.wavelengthColor, DisplayName="pseudo-wavelength (km)");
    clabel(C,h,options.wavelengths,'Color',options.wavelengthColor,'LabelSpacing',options.labelSpacing)
    hold off
    % y axis ticks
    % select reasonable jMode spacing for ticks
    jMode = [0:9,10:10:90,100:100:900];
    jMode = jMode(jMode<=max(wvt.j)); % truncate at maximum mode
    jModePseudo = jMode;
    jModePseudo(1) = jModePseudo(2)/2;
    % set tick locations
    set(ax, 'YTick', jModePseudo);
    set(ax, 'YMinorTick','off')
    % set tick labels
    maxLeadDigit = 5;
    leadDigit = floor(jMode ./ 10.^floor(log10(jMode)));
    leadDigit(1) = 0;
    labels = strings(size(jMode));
    labels(leadDigit <= maxLeadDigit) = string(jMode(leadDigit <= maxLeadDigit));
    set(ax,'YTickLabel', labels);

    % wave energy spectrum
    ax = nexttile(tl_inner,2);
    val = log10(TE_Apm_j_kR);
    pcolor(ax,2*pi./kPseudoLocation/1000,jPseudoLocation,val), shading flat,
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    axis square
    cb = colorbar('southoutside');
    clim(ax,energy_limits);
    title('wave energy spectrum')
    xlabel('horizontal wavelength (km)')
    ylabel("Vertical mode number")
    set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
    colormap(ax, cmLinMono);
    cb.Label.String = "log10(m^3 s^{-2})";
    set(gca,'layer','top'),
    hold on
    % add frequency contours. flipud/fliplr gives nicer clabel placement.
    % Note: have to pad quantities to work right with our log-log axes. And here
    % have to "double pad" to work right with pcolor's box shift and
    % reversed axis direction. 
    omegaPadded = cat(1,[self.omega_jk(1,:);self.omega_jk(1,:)],self.omega_jk(1:end-1,:));
    omegaPadded = cat(2,[omegaPadded(:,1) omegaPadded(:,1)],omegaPadded(:,1:end-1));
    omegaJK = interpn(KPaddedWVD,JPaddedWVD,omegaPadded.',KPseudoLocation,JPseudoLocation,"linear");
    [C,h] = contour(ax,flipud(2*pi./kPseudoLocation/1000),jPseudoLocation,fliplr(omegaJK'/wvt.f),options.frequencies,'LineWidth',options.lineWidth,'Color',options.frequencyColor, DisplayName="frequency (f)", HandleVisibility='off');
    clabel(C,h,options.frequencies,'Color',options.frequencyColor,'LabelSpacing',options.labelSpacing)
    % add pseudoWavelength. flipud/fliplr gives nicer clabel placement.
    [C,h] = contour(ax,flipud(2*pi./kPseudoLocation/1000),jPseudoLocation,fliplr(2*pi./KPseudoRadial'/1000),options.wavelengths,'LineWidth',options.lineWidth,'Color',options.wavelengthColor, DisplayName="pseudo-wavelength (km)");
    clabel(C,h,options.wavelengths,'Color',options.wavelengthColor,'LabelSpacing',options.labelSpacing)
    hold off
    % y axis ticks
    % select reasonable jMode spacing for ticks
    jMode = [0:9,10:10:90,100:100:900];
    jMode = jMode(jMode<=max(wvt.j)); % truncate at maximum mode
    jModePseudo = jMode;
    jModePseudo(1) = jModePseudo(2)/2;
    % set tick locations
    set(ax, 'YTick', jModePseudo);
    set(ax, 'YMinorTick','off')
    % set tick labels
    maxLeadDigit = 5;
    leadDigit = floor(jMode ./ 10.^floor(log10(jMode)));
    leadDigit(1) = 0;
    labels = strings(size(jMode));
    labels(leadDigit <= maxLeadDigit) = string(jMode(leadDigit <= maxLeadDigit));
    set(ax,'YTickLabel', labels);

end

% Add a large label to the left of the first row (row 1)
% textAxes = axes('Position', [0, 0, 1, 1], 'Visible', 'off'); % Dummy invisible axes
% text(0.07, 0.74, 'Geostrophic', 'FontSize', 12, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label
% text(0.07, 0.34, 'Wave', 'FontSize', 12, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label

end