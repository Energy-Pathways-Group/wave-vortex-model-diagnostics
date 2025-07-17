function fig = plotFluidDecompositionMultipanel(self,options)
arguments
    self WVDiagnostics
    options.visible = "on"
    options.iTime
    options.title
    options.yForXZSlice
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

if ~isfield(options,"title")
    options.title = sprintf('%d days',round(wvt.t/86400));
end

cmDivRWB = self.cmocean('balance'); % diverging positive-negative

% set limits
zeta_limits = [-0.2 0.2];

% location for x-z section
if ~isfield(options,"yForXZSlice")
    iY = round(wvt.Nx/2);
else
    iY = round(options.yForXZSlice/(wvt.y(2)-wvt.y(1)));
end

nColumns = 3;
if nColumns == 3
    figPos = [50 50 900 600];
else
    figPos = [50 50 600 615];
end

fig = figure(Units='points',Position=figPos,Visible = options.visible);
set(gcf,'PaperPositionMode','auto')


tl = tiledlayout(2,nColumns,TileSpacing="tight");

if options.title ~= "none"
    title(tl, options.title, 'Interpreter', 'none')
end

% wave and geosgrophic vorticity v_x - u_y
zeta_z_g = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);  % geostrophic
zeta_z_w = wvt.diffX(wvt.v_w) - wvt.diffY(wvt.u_w);  % wave

    function makeVorticityXZPlot(zeta_z)
        val = squeeze(zeta_z(:,iY,:)/wvt.f);
        pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
        axis square
        colormap(ax, cmDivRWB);
        set(gca,'Layer','top','TickLength',[0.015 0.015])
        clim(ax, zeta_limits);
    end

    function makeVorticityXYPlot(zeta_z)
        val = zeta_z(:,:,end)/wvt.f;
        pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
        hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,Color=0*[1 1 1],LineWidth=2); % add line for x-z section
        axis square
        colormap(ax, cmDivRWB);
        set(gca,'Layer','top','TickLength',[0.015 0.015])
        clim(ax, zeta_limits);
    end

% geostrophic vorticity section
ax = nexttile(tl,1);
makeVorticityXYPlot(wvt.zeta_z);
xticklabels([])
ylabel('y-distance (km)')
title(ax, "total")

ax = nexttile(tl,2);
makeVorticityXYPlot(zeta_z_w)
xticklabels([])
yticklabels([])
title(ax, "wave")

ax = nexttile(tl,3);
makeVorticityXYPlot(zeta_z_g);
xticklabels([])
yticklabels([])
title(ax, "geostrophic")
cb = colorbar("eastoutside");
cb.Label.String = "$\zeta$/f";
cb.Label.Interpreter = 'latex';

ax = nexttile(tl,4);
makeVorticityXZPlot(wvt.zeta_z);
xlabel('x-distance (km)')
ylabel('depth (km)')

ax = nexttile(tl,5);
makeVorticityXZPlot(zeta_z_w)
xlabel('x-distance (km)')
yticklabels([])

ax = nexttile(tl,6);
makeVorticityXZPlot(zeta_z_g);
xlabel('x-distance (km)')
yticklabels([])
cb = colorbar("eastoutside");
cb.Label.String = "$\zeta$/f";
cb.Label.Interpreter = 'latex';

end