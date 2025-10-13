function [sources, sinks, inertial, ddt, energy] = filterEnergyForSourcesSinksReservoirs2(self,options)
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

    kp = self.kPseudoRadial;
    forcing_fluxes_kp = self.quadraticEnergyFluxesTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);
    forcing_exact_kp = self.exactEnergyFluxesTemporalAverage(timeIndices=options.timeIndices);

    % First convert the forcing fluxes to kPseudoRadial
    for iForce=1:length(forcing_fluxes_kp)
        forcing_fluxes_kp(iForce).te_gmda = self.transformToPseudoRadialWavenumberA0(forcing_fluxes_kp(iForce).te_gmda);
        forcing_fluxes_kp(iForce).te_wave = self.transformToPseudoRadialWavenumberApm(forcing_fluxes_kp(iForce).te_wave);
    end
    for iForce=1:length(forcing_exact_kp)
        forcing_exact_kp(iForce).te = self.transformToPseudoRadialWavenumberA0(forcing_exact_kp(iForce).te);
    end

    % Now find which indices are outside the damping region
    damping_index = [forcing_fluxes_kp.name] == "adaptive_damping";
    if ~any(damping_index)
        error("Unable to find adaptive damping and thus do not know how to proceed. I guess we could look for another closure.");
    end
    damping_kp = forcing_fluxes_kp(damping_index).te_gmda + forcing_fluxes_kp(damping_index).te_wave;
    NoDamp = abs(cumsum(damping_kp/self.flux_scale)) < options.fluxTolerance;

    % Sort the forcing (both quadratic and exact) into damped and undamped
    forcing = self.filterFluxesForReservoir(forcing_fluxes_kp,filter=@(v) sum(sum(v(NoDamp))));
    forcing_fluxes_damp = self.filterFluxesForReservoir(forcing_fluxes_kp,filter=@(v) sum(sum(v(~NoDamp))));
    forcing_exact = self.filterFluxesForReservoir(forcing_exact_kp,filter=@(v) sum(sum(v(NoDamp))));
    forcing_exact_damp = self.filterFluxesForReservoir(forcing_exact_kp,filter=@(v) sum(sum(v(~NoDamp))));

    % create a new reservoir, te_damp, which contains the combined geo and
    % wave forcing in the damped region. Same for the exact flux.
    for iForce=1:length(forcing)
        forcing(iForce).te_damp = forcing_fluxes_damp(iForce).te_gmda + forcing_fluxes_damp(iForce).te_wave;
        forcing(iForce).te_exact = forcing_exact(iForce).te;
        forcing(iForce).te_exact_damp = forcing_exact_damp(iForce).te;
    end

    % inertial_triads = self.filterFluxesForReservoir(inertial_jk,filter=@(v) sum(sum(v(NoDamp))));
    
    % Now deal with the inertial triads. Here we use the *sparse* pseudo
    % radial axis, so we need to find the appropriate indices again.
    [inertial_fluxes_g_kps, inertial_fluxes_w_kps, kps] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(timeIndices=options.timeIndices);
    NoDampKps = kps <= max(kp(NoDamp));

    gmda_tx_wave = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-wwg").flux(NoDampKps) + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-ggw").flux(NoDampKps));
    wave_tx_gmda = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-wwg").flux(NoDampKps) + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-ggw").flux(NoDampKps));
    gmda_tx_wave = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-wwg").flux + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-ggw").flux);
    wave_tx_gmda = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-wwg").flux + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-ggw").flux);

    % There is one route for te_gmda to transfer to te_wave, and that is
    % the direct transfer terms summed. Although, *some* of that energy
    % might land in the damping region of the waves, and thus actually land
    % in te_damp. The two cascade terms will also transfer to te_damp.
    inertial(1).name = "te_gmda";
    inertial(1).fancyName = WVDiagnostics.fancyNameForName(inertial(1).name);
    inertial(1).te_gmda = 0;
    inertial(1).te_damp = inertial_triads(1).te_gmda + inertial_triads(2).te_gmda + inertial_triads(3).te_gmda + inertial_triads(1).te_wave;

    inertial(2).name = "te_wave";
    inertial(2).fancyName = WVDiagnostics.fancyNameForName(inertial(2).name);
    inertial(2).te_damp = inertial_triads(4).te_wave + inertial_triads(2).te_wave + inertial_triads(3).te_wave + inertial_triads(4).te_gmda;
    inertial(2).te_wave = 0;

    inertial(3).name = "te_exact";
    inertial(3).fancyName = "exact undamped reservoir";
    inertial(3).te_gmda = 0;
    inertial(3).te_damp = forcing(1).te_exact;
    inertial(3).te_wave = 0;

    % remove nonlinear advection, now that we copied the values we needed
    forcing(1) = [];

    % divide into sources and sinks
    iSink = 1; iSource = 1;
    for iForce=1:length(forcing)
        if forcing(iForce).te_exact + forcing(iForce).te_exact_damp < 0
            sinks(iSink) = forcing(iForce); iSink = iSink + 1;
        else
            sources(iSource) = forcing(iForce); iSource = iSource + 1;
        end
    end

    % Alt 1
    inertial(1).te_wave = inertial_triads(4).te_gmda - inertial_triads(1).te_wave;
    inertial(2).te_gmda = inertial_triads(1).te_wave - inertial_triads(4).te_gmda;

    [reservoirEnergy, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices,shouldIncludeExactTotalEnergy=true);
    

    % ape_bg = self.diagfile.readVariables("ape_bg");
    % ape_bg = ape_bg(options.timeIndices);
    % ddt.te_bg = (ape_bg(end) - ape_bg(1))/(t(end)-t(1));
    energy.te_gmda = mean(reservoirEnergy(1).energy);
    energy.te_wave = mean(reservoirEnergy(2).energy);
    energy.te_exact = mean(reservoirEnergy(3).energy);
    E = self.wvt.transformToRadialWavenumber(self.wvt.Apm_TE_factor.*(abs(self.wvt.Ap).^2 + abs(self.wvt.Am).^2) + self.wvt.A0_TE_factor.*(abs(self.wvt.A0).^2));
    E_big = zeros(length(self.j),length(self.kRadial));
    E_big(1:size(E,1),1:size(E,2),:) = E;
    energy.te_damp = sum( E_big(:) .* ~NoDamp(:));
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