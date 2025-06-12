function fig = plotFluidStateMultipanel(self,options)
arguments
    self WVDiagnostics
    options.visible = "on"
    options.iTime
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

cmDivRWB = self.cmocean('balance'); % diverging positive-negative
cmLinMono = self.cmocean('dense') ; % linear monochrome for data>0

% set limits
zeta_limits = [-0.2 0.2];
energy_limits = [-8 0];
enstrophy_limits = [-16 -9];

% create some nice tick labels to show wavelength
% ticks_x = [500;300;200;100;50;30;20;10];
% ticks_x = [500;50;30;20;15;12;10];
ticks_x = round(2*pi./(1000.*linspace(wvt.kRadial(2),wvt.kRadial(end),6)));
labels_x = cell(length(ticks_x),1);
for i=1:length(ticks_x)
    labels_x{i} = sprintf('%.0f',ticks_x(i));
end
ticks_x = 2*pi./(1e3*ticks_x);

% color for some annotations
textColor = '[.5,.5,.5]';

% location for x-z section
iY = round(wvt.Nx/2);

% create the lines of constant frequency
[omegaN,n] = wvt.transformToRadialWavenumber(abs(wvt.Omega),ones(size(wvt.Omega)));
omegaJK = (omegaN./n)/wvt.f;

% create the lines of constant deformation radius
deformationJK = repmat(sqrt(wvt.Lr2)./1000,1,length(wvt.kRadial));

fig = figure(Units='points',Position=[50 50 900 600],Visible = options.visible);
set(gcf,'PaperPositionMode','auto')

tl = tiledlayout(2,3,TileSpacing="tight");

title(tl, sprintf(':  %d days',round(wvt.t/86400)), 'Interpreter', 'none')

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

% geostrophic vorticity section
ax = nexttile(tl,1);
val = squeeze(zeta_z_g(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
% title("geostrophic x-z vorticity")
title("x-z vorticity")
axis square
colormap(ax, cmDivRWB);
xticklabels([])
ylabel('Depth (km)')
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);

% wave vorticity section
ax = nexttile(tl,4);
val = squeeze(zeta_z_w(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
% title("wave x-z vorticity")
% title("x-z vorticity")
axis square
colormap(ax, cmDivRWB);
xlabel('x-distance (km)')
ylabel('Depth (km)')
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);
cb = colorbar("southoutside");
cb.Label.String = "$\zeta$/f";
cb.Label.Interpreter = 'latex';

% geostrophic surface vorticity
ax = nexttile(tl,2);
val = zeta_z_g(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:'); % add line for x-z section
title("surface vorticity")
% title("geostrophic surface vorticity")
axis square
xticklabels([])
yticklabels([])
% ylabel('y-distance (km)')
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
colormap(ax, cmDivRWB);
clim(ax, zeta_limits);

% wave surface vorticity
ax = nexttile(tl,5);
val = zeta_z_w(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:'); % add line for x-z section
% title("surface vorticity")
% title("wave surface vorticity")
axis square
colormap(ax, cmDivRWB);
xlabel('x-distance (km)')
yticklabels([])
% ylabel('y-distance (km)')
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);
cb = colorbar("southoutside");
cb.Label.String = "$\zeta$/f";
cb.Label.Interpreter = 'latex';

% geostrophic energy spectrum
ax = nexttile(tl,3);
val = log10(TE_A0_j_kR);
pcolor(ax,wvt.kRadial,wvt.j,val), shading flat,
% title('geostrophic energy')
title('energy spectrum')
axis square
clim(ax,energy_limits);
xticks(ticks_x)
xticklabels([])
% xticklabels(labels_x)
% xlabel('wavelength (km)')
% xlabel("wavenumber")
ylabel("Vertical mode")
set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
colormap(ax, cmLinMono);
% cb.Label.String = "log10(m^3 s^{-2})";
% create some nice tick labels to show deformation radius
yticksTemp = yticks;
ticks_y = sqrt(wvt.Lr2)./1000;
labels_y = cell(length(yticksTemp),1);
for i=1:length(yticksTemp)
    labels_y{i} = sprintf('%0.1f',ticks_y(yticksTemp(i)+1));
end
text(1.2*max(xlim)*ones(size(yticksTemp)),yticksTemp,labels_y,'Color',textColor,'HorizontalAlignment','left')
text(1.3*max(xlim),1.05*max(ylim),'L_r (km)','Color',textColor,'HorizontalAlignment','right')
% % add deformation radius contours
% set(gca,'layer','top'),
% hold on
% v = round(sqrt(wvt.Lr2(1:2:end))./1000*10)/10;
% [C,h] = contour(wvt.kRadial,wvt.j',deformationJK,v,'LineWidth',1,'Color',textColor);
% clabel(C,h,v,'Color',textColor,'LabelSpacing',400)

% wave energy spectrum
ax = nexttile(tl,6);
val = log10(TE_Apm_j_kR);
pcolor(ax,wvt.kRadial,wvt.j,val), shading flat,
% title('wave energy')
% title('energy spectrum')
axis square
cb = colorbar('southoutside');
clim(ax,energy_limits);
xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)')
clim(ax,energy_limits);
% xlabel("wavenumber (rad m^{-1})")
% xticklabels([])
ylabel("Vertical mode")
set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
colormap(ax, cmLinMono);
cb.Label.String = "log10(m^3 s^{-2})";
% add frequency contours
set(gca,'layer','top'),
hold on
v = [1.01 1.05 1.2 1.5 2 4 8 16];
[C,h] = contour(wvt.kRadial(2:end),wvt.j',(omegaJK(:,2:end)),v,'LineWidth',1,'Color',textColor);
clabel(C,h,v,'Color',textColor,'LabelSpacing',400)

% Add a large label to the left of the first row (row 1)
textAxes = axes('Position', [0, 0, 1, 1], 'Visible', 'off'); % Dummy invisible axes
text(0.07, 0.74, 'Geostrophic', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label
text(0.07, 0.34, 'Wave', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label

end