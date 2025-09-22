function [forcing, inertial, ddt] = filterEnergyForSourcesSinksReservoirs(self,options)
% This function returns values assuming three reservoirs: geo, wave, and
% damping. The damping resevoir is just scales below a threshold, wave or
% geostrophic. It also returns the exact and exact-damp resevoirs. 
%
% The forcing struct has the the forcing on the three/two different
% reservoirs
%
% The inertial struct has the flux from the two reservoirs (wave
% geostrophic) to each other and to the damping region.
%
% The forcing struct also include the nonlinear advection, which has the
% flux to the damping region.
%
% The ddt struct contains the change in total energy, closing the energy
% budget.
arguments
    self WVDiagnostics
    options.customNames = configureDictionary("string","string")
    options.fluxTolerance = 1e-2;
    options.timeIndices = Inf;
    options.shouldShowReservoirEnergy = true
    options.shouldShowExactValues = true
    options.shouldSeparateClosureRegion = true
end
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave];

    forcing_fluxes_jk = self.quadraticEnergyFluxesTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);
    forcing_exact_jk = self.exactEnergyFluxesTemporalAverage(timeIndices=options.timeIndices);
    inertial_jk = self.quadraticEnergyTriadFluxesTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

    % We could be smarter about this. But roughly speaking we want a region
    % where none of the forcing shows up above our precision reporting
    % level. It turns out our tidal forcing shows up if I do no_damp, so
    % gotta put in the middle to strike the balance.
    % Here were some failed ideas... the current idea is to split the
    % difference, literally average the two values. It seems to work well.
    
    svv = self.wvt.forcingWithName("adaptive damping");
    [J,K] = ndgrid(self.j,self.kRadial);

    % NoDamp = (K < svv.k_damp) & (J < svv.j_damp);
    % NoDamp = (K < svv.k_no_damp) & (J < svv.j_no_damp);
    NoDamp = (K < (svv.k_damp+svv.k_no_damp)/2) & (J < (svv.j_damp + svv.j_no_damp)/2);
    forcing = self.filterFluxesForReservoir(forcing_fluxes_jk,filter=@(v) sum(sum(v(NoDamp))));
    forcing_fluxes_damp = self.filterFluxesForReservoir(forcing_fluxes_jk,filter=@(v) sum(sum(v(~NoDamp))));
    forcing_exact = self.filterFluxesForReservoir(forcing_exact_jk,filter=@(v) sum(sum(v(NoDamp))));
    forcing_exact_damp = self.filterFluxesForReservoir(forcing_exact_jk,filter=@(v) sum(sum(v(~NoDamp))));
    inertial_triads = self.filterFluxesForReservoir(inertial_jk,filter=@(v) sum(sum(v(NoDamp))));

    for iForce=1:length(forcing)
        forcing(iForce).te_damp = forcing_fluxes_damp(iForce).te_gmda + forcing_fluxes_damp(iForce).te_wave;
        % forcing(iForce).te_quadratic = forcing(iForce).te_gmda + forcing(iForce).te_wave + forcing(iForce).te_damp;
        forcing(iForce).te_exact = forcing_exact(iForce).te;
        forcing(iForce).te_exact_damp = forcing_exact_damp(iForce).te;
    end

    inertial(1).name = "te_gmda";
    inertial(1).fancyName = WVDiagnostics.fancyNameForName(inertial(1).name);
    inertial(1).te_gmda = 0;
    inertial(1).te_damp = inertial_triads(1).te_gmda + inertial_triads(2).te_gmda + inertial_triads(3).te_gmda + inertial_triads(1).te_wave;

    inertial(2).name = "te_wave";
    inertial(2).fancyName = WVDiagnostics.fancyNameForName(inertial(2).name);
    inertial(2).te_damp = inertial_triads(4).te_wave + inertial_triads(2).te_wave + inertial_triads(3).te_wave + inertial_triads(4).te_gmda;
    inertial(2).te_wave = 0;

    % Alt 1
    inertial(1).te_wave = inertial_triads(4).te_gmda - inertial_triads(1).te_wave;
    inertial(2).te_gmda = inertial_triads(1).te_wave - inertial_triads(4).te_gmda;

    [reservoirEnergy, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices,shouldIncludeExactTotalEnergy=true);
    

    % ape_bg = self.diagfile.readVariables("ape_bg");
    % ape_bg = ape_bg(options.timeIndices);
    % ddt.te_bg = (ape_bg(end) - ape_bg(1))/(t(end)-t(1));
    ddt.te_gmda = (reservoirEnergy(1).energy(end) - reservoirEnergy(1).energy(1))/(t(end)-t(1));
    ddt.te_wave = (reservoirEnergy(2).energy(end) - reservoirEnergy(2).energy(1))/(t(end)-t(1));
    ddt.te_exact = (reservoirEnergy(3).energy(end) - reservoirEnergy(3).energy(1))/(t(end)-t(1));

end

% One can replace "inertial_triads(1).te_wave" with E0_ggw
% and replace "inertial_triads(4).te_gmda" with Epm_wwg
% [E0_ggw_jkt,Epm_wwg_jkt] = self.quadraticEnergyMirrorTriadsUndamped(timeIndices=options.timeIndices);
% E0_ggw = sum(sum(mean(NoDamp.*E0_ggw_jkt,3),2),1);
% Epm_wwg = sum(sum(mean(NoDamp.*Epm_wwg_jkt,3),2),1);
% E0_ggw = sum(sum(mean(E0_ggw_jkt,3),2),1);
% Epm_wwg = sum(sum(mean(Epm_wwg_jkt,3),2),1);
% fprintf("low-pass ggw: %f, total ggw: %f\n", E0_ggw/self.flux_scale,inertial_triads(1).te_wave/self.flux_scale);
% fprintf("low-pass wwg: %f, total wwg: %f\n", Epm_wwg/self.flux_scale,inertial_triads(4).te_gmda/self.flux_scale);
% some geostrophic energy goes straight to the damped wave modes
% % Alt 2: asymmetric. Not good
% inertial(1).te_wave = inertial_triads(4).te_gmda - E0_ggw;
% inertial(2).te_gmda = inertial_triads(1).te_wave - Epm_wwg;
%
% % Alt 3:
% inertial(1).te_wave = Epm_wwg - E0_ggw;
% inertial(2).te_gmda = E0_ggw - Epm_wwg;