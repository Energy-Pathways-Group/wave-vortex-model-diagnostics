function fig = plotEnergyFlux1D(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
    options.visible = "on"
end



[inertial_fluxes_g, inertial_fluxes_w, kp] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(timeIndices=options.timeIndices);


radialWavelength = 2*pi./kp/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

filter = @(v) cumsum(v)/self.flux_scale;

fig = figure(Visible=options.visible);
tl = tiledlayout(fig,GridSize=[2 1],TileSpacing='tight');

ax = nexttile(tl);
plot(radialWavelength,filter(inertial_fluxes_g(1).te_gmda)), xscale('log'), hold on
for i=2:length(inertial_fluxes_g)
    plot(radialWavelength,filter(inertial_fluxes_g(i).te_gmda))
end
set(ax,'XDir','reverse')

ax = nexttile(tl);
plot(radialWavelength,filter(inertial_fluxes_w(1).te_wave)), xscale('log'), hold on
for i=2:length(inertial_fluxes_w)
    plot(radialWavelength,filter(inertial_fluxes_w(i).te_wave))
end
set(ax,'XDir','reverse')

end