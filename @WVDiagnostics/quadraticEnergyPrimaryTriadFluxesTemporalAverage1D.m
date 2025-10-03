function [inertial_fluxes_g, inertial_fluxes_w, kp] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

[M_wwg, M_ggw, kp] = self.quadraticEnergyMirrorTriadFluxes1D(timeIndices=options.timeIndices);

flux_interp = @(v) cat(1,zeros(1,size(v,2)),diff(interp1(self.kPseudoRadial,cumsum(v),kp)));

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_g = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.geostrophic_mda,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_g)
    val = self.transformToPseudoRadialWavenumber(EnergyReservoir.geostrophic_mda,fluxes_g(idx).te_gmda);
    fluxes_g(idx).flux = flux_interp(val);
end

fluxes_w = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.wave,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_w)
    val = self.transformToPseudoRadialWavenumber(EnergyReservoir.wave,fluxes_w(idx).te_wave);
    fluxes_w(idx).flux = flux_interp(val);
end

inertial_fluxes_g(1).flux = fluxes_g([fluxes_g.name] == "gmda_gmda").flux;
inertial_fluxes_g(1).name = "ggg";
inertial_fluxes_g(1).fancyName = "$[g{\nabla}g]_g$";

inertial_fluxes_g(2).flux = fluxes_g([fluxes_g.name] == "gmda_wave").flux + fluxes_g([fluxes_g.name] == "wave_gmda").flux + M_ggw;
inertial_fluxes_g(2).name = "ggw";
inertial_fluxes_g(2).fancyName = "$[g{\nabla}w]_g+[w{\nabla}g]_g+\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_g(3).flux = fluxes_g([fluxes_g.name] == "wave_wave").flux;
inertial_fluxes_g(3).name = "tx-wwg";
inertial_fluxes_g(3).fancyName = "$[w{\nabla}w]_g$";

inertial_fluxes_g(4).flux = -M_ggw;
inertial_fluxes_g(4).name = "tx-ggw";
inertial_fluxes_g(4).fancyName = "$-\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_w(1).flux = fluxes_w([fluxes_w.name] == "wave_wave").flux;
inertial_fluxes_w(1).name = "www";
inertial_fluxes_w(1).fancyName = "$[w{\nabla}w]_w$";

inertial_fluxes_w(2).flux = fluxes_w([fluxes_w.name] == "gmda_wave").flux + fluxes_w([fluxes_w.name] == "wave_gmda").flux + M_wwg;
inertial_fluxes_w(2).name = "wwg";
inertial_fluxes_w(2).fancyName = "$[g{\nabla}w]_w+[w{\nabla}g]_w+\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(3).flux = -M_wwg;
inertial_fluxes_w(3).name = "tx-wwg";
inertial_fluxes_w(3).fancyName = "$-\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(4).flux = fluxes_w([fluxes_w.name] == "gmda_gmda").flux;
inertial_fluxes_w(4).name = "tx-ggw";
inertial_fluxes_w(4).fancyName = "$[g{\nabla}g]_w$";

end