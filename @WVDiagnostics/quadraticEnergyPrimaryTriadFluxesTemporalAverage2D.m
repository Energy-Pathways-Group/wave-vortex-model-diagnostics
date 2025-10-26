function [inertial_fluxes_g, inertial_fluxes_w, ks, js] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

[J,K] = ndgrid(self.jWavenumber,self.kRadial);
js = self.sparseJWavenumberAxis;
ks = self.sparseKRadialAxis;
[Js,Ks] = ndgrid(js,ks);
matrixSize = [length(js) length(ks)];

flux_interp = @(v) diff(diff( cat(2,zeros(length(js)+1,1),cat(1,zeros(1,length(ks)),interpn(J,K,cumsum(cumsum(v,1),2),Js,Ks))), 1,1 ),1,2);

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_g = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.geostrophic_mda,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_g)
    fluxes_g(idx).flux = flux_interp(fluxes_g(idx).te_gmda);
end

fluxes_w = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.wave,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_w)
    fluxes_w(idx).flux = flux_interp(fluxes_w(idx).te_wave);
end

if ~self.diagfile.hasGroupWithName("mirror-flux-2d-wwg")
    M_wwg = zeros(matrixSize);
    fprintf("Did not find the 2D mirror fluxes for wwg, assuming it is zero.\n");
else
    M_wwg = mean(self.quadraticEnergyMirrorTriadFluxes2D(timeIndices=options.timeIndices,mirrorTriad="wwg"),3);
end

if ~self.diagfile.hasGroupWithName("mirror-flux-2d-ggw")
    M_ggw = zeros(matrixSize);
    fprintf("Did not find the 2D mirror fluxes for ggw, assuming it is zero.\n");
else
    M_ggw = mean(self.quadraticEnergyMirrorTriadFluxes2D(timeIndices=options.timeIndices,mirrorTriad="ggw"),3);
end

inertial_fluxes_g(1).flux = fluxes_g([fluxes_g.name] == "gmda_gmda").flux;
inertial_fluxes_g(1).name = "ggg";
inertial_fluxes_g(1).fancyName = "ggg cascade";
% inertial_fluxes_g(1).fancyName = "$[g{\nabla}g]_g$";


inertial_fluxes_g(2).flux = fluxes_g([fluxes_g.name] == "gmda_wave").flux + fluxes_g([fluxes_g.name] == "wave_gmda").flux + M_ggw;
inertial_fluxes_g(2).name = "ggw";
inertial_fluxes_g(2).fancyName = "ggw cascade";
% inertial_fluxes_g(2).fancyName = "$[g{\nabla}w]_g+[w{\nabla}g]_g+\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_g(3).flux = fluxes_g([fluxes_g.name] == "wave_wave").flux;
inertial_fluxes_g(3).name = "tx-wwg";
inertial_fluxes_g(3).fancyName = "wwg transfer";
% inertial_fluxes_g(3).fancyName = "$[w{\nabla}w]_g$";

inertial_fluxes_g(4).flux = -M_ggw;
inertial_fluxes_g(4).name = "tx-ggw";
inertial_fluxes_g(4).fancyName = "ggw transfer";
% inertial_fluxes_g(4).fancyName = "$-\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_w(1).flux = fluxes_w([fluxes_w.name] == "wave_wave").flux;
inertial_fluxes_w(1).name = "www";
inertial_fluxes_w(1).fancyName = "www cascade";
% inertial_fluxes_w(1).fancyName = "$[w{\nabla}w]_w$";

inertial_fluxes_w(2).flux = fluxes_w([fluxes_w.name] == "gmda_wave").flux + fluxes_w([fluxes_w.name] == "wave_gmda").flux + M_wwg;
inertial_fluxes_w(2).name = "wwg";
inertial_fluxes_w(2).fancyName = "wwg cascade";
% inertial_fluxes_w(2).fancyName = "$[g{\nabla}w]_w+[w{\nabla}g]_w+\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(3).flux = -M_wwg;
inertial_fluxes_w(3).name = "tx-wwg";
inertial_fluxes_w(3).fancyName = "wwg transfer";
% inertial_fluxes_w(3).fancyName = "$-\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(4).flux = fluxes_w([fluxes_w.name] == "gmda_gmda").flux;
inertial_fluxes_w(4).name = "tx-ggw";
inertial_fluxes_w(4).fancyName = "ggw transfer";
% inertial_fluxes_w(4).fancyName = "$[g{\nabla}g]_w$";

end