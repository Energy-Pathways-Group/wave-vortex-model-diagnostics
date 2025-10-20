function addTriadFluxesForReservoirGroupAtTime(self,options)
arguments
    self WVDiagnostics
    options.triadVar
    options.flowComponents
    options.wvt
    options.outputIndex
end
triadVar = options.triadVar;
flowComponents = options.flowComponents;
wvt = options.wvt;
outputIndex = options.outputIndex;

    for i=1:length(flowComponents)
        for j=1:i
            if i==j
                [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(flowComponents(i),flowComponents(j));
            else
                [Fp_ij,Fm_ij,F0_ij] = wvt.nonlinearFluxForFlowComponents(flowComponents(i),flowComponents(j));
                [Fp_ji,Fm_ji,F0_ji] = wvt.nonlinearFluxForFlowComponents(flowComponents(j),flowComponents(i));
                Fp = Fp_ji + Fp_ij;
                Fm = Fm_ji + Fm_ij;
                F0 = F0_ji + F0_ij;
            end
            [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
            for k=1:length(flowComponents)
                var = triadVar{"T_" + i + "_" + j + "_" + k};
                dE = sum(flowComponents(k).maskAp(:).*Ep(:) + flowComponents(k).maskAm(:).*Em(:) + flowComponents(k).maskA0(:).*E0(:));
                var.setValueAlongDimensionAtIndex(dE,'t',outputIndex);
            end
        end
    end

end