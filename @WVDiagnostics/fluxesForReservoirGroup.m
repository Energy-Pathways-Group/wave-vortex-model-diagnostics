function [transferFlux, forcingFlux, ddt, energy] = fluxesForReservoirGroup(self,options)
arguments
    self WVDiagnostics
    options.outputfile NetCDFGroup
    options.name string = "reservoir-damped-wave-geo"
    options.timeIndices
    options.outputFormat string {mustBeMember(options.outputFormat,{'matrix','struct','table'})} = 'matrix'
end

if ~isfield(options,"outputfile")
    options.outputfile = self.diagfile;
end

t = self.diagfile.readVariables("t");
if ~isfield(options,"timeIndices")
    timeIndices = 1:length(t);
else
    timeIndices = options.timeIndices;
end

if ~options.outputfile.hasGroupWithName(options.name)
    error("Unable to find the group "+options.name+".");
end
group = options.outputfile.groupWithName(options.name);
flowComponentNames = reshape(group.attributes('flow-components'),[],1);
flowComponentNamesNormalized = replace(replace(flowComponentNames,"-","_")," ","_");

forcingNames = self.forcingNames;
forcingNamesNormalized = replace(replace(self.forcingNames,"-","_")," ","_");

forcingFlux = zeros(length(forcingNames),length(flowComponentNames));
for iForce=1:length(forcingNames)
    for k=1:length(flowComponentNames)
        dE = group.readVariables(forcingNamesNormalized(iForce)+"_"+k);
        forcingFlux(iForce,k) = mean(dE(timeIndices));
    end
end
% % remove nonlinear advection
% forcingNames(1) = [];
% forcingNamesNormalized(1) = [];
% forcingFlux(1,:) = [];

ddt = zeros(1,length(flowComponentNames));
energy = zeros(1,length(flowComponentNames));
for k=1:length(flowComponentNames)
    E = group.readVariables("E_"+k);
    energy(k) = mean(E(timeIndices));
    ddt(k) = (E(timeIndices(end)) - E(timeIndices(1)))/(t(timeIndices(end)) - t(timeIndices(1)));
end

% loop over T_i_j_k
transferFlux = zeros(length(flowComponentNames),length(flowComponentNames));
for k=1:length(flowComponentNames)
    % Flux *into* the k-th reservoir
    for i=1:length(flowComponentNames)
        for j=1:i
            
            if i==j
                % has the form T_i_i_k OR T_k_k_k
                E = group.readVariables("T_" + i + "_" + j + "_" + k);
                E = mean(E(timeIndices));
                transferFlux(i,k) = transferFlux(i,k) + E;
            elseif i==k % our for-loop is such that i will reach k, but never let j get that high
                % has the form T_k_j_k, so we need -T_k_k_j
                E = group.readVariables("T_" + k + "_" + k + "_" + j);
                E = mean(E(timeIndices));
                transferFlux(j,k) = transferFlux(j,k) - E;
            elseif j==k % our for-loop is such that i will reach k, but never let j get that high
                % has the form T_i_k_k, so we need -T_k_k_i
                E = group.readVariables("T_" + k + "_" + k + "_" + i);
                E = mean(E(timeIndices));
                transferFlux(i,k) = transferFlux(i,k) - E;
            else
                % now we have a true mixed triad
                E_i_j_k = group.readVariables("T_" + i + "_" + j + "_" + k);
                E_j_k_i = group.readVariables("T_" + max(j,k) + "_" + min(j,k) + "_" + i);
                E_k_i_j = group.readVariables("T_" + max(i,k) + "_" + min(i,k) + "_" + j);

                E_i_j_k = mean(E_i_j_k(timeIndices));
                E_j_k_i = mean(E_j_k_i(timeIndices));
                E_k_i_j = mean(E_k_i_j(timeIndices));

                % d = (E_j_k_i + E_k_i_j + E_i_j_k)/self.flux_scale; % / max([abs(E_j_k_i),abs(E_k_i_j),abs(E_i_j_k)]);
                % disp( "T_" + i + "_" + j + "_" + k + " with relative error " + string(d))

                % a -> i, b -> j, c -> k
                [Tji, Tki, Tkj] = transfersFromFluxes(E_j_k_i,E_k_i_j,E_i_j_k);
                transferFlux(i,k) = transferFlux(i,k) - Tki;
                transferFlux(j,k) = transferFlux(j,k) - Tkj;
            end
        end
    end
end

if options.outputFormat == "struct" || options.outputFormat == "table"
    forcingFluxStruct(length(forcingNames)) = struct("name","placeholder");
    for iForce=1:length(forcingNames)
        forcingFluxStruct(iForce).name = forcingNamesNormalized(iForce);
        forcingFluxStruct(iForce).fancyName = forcingNames(iForce);
        for k=1:length(flowComponentNames)
            forcingFluxStruct(iForce).(flowComponentNamesNormalized(k)) = forcingFlux(iForce,k);
        end
    end
    forcingFlux = forcingFluxStruct;

    transferFluxStruct(length(flowComponentNames)) = struct("name","placeholder");
    for i=1:length(flowComponentNames)
        transferFluxStruct(i).name = flowComponentNamesNormalized(i);
        transferFluxStruct(i).fancyName = flowComponentNames(i);
        for k=1:length(flowComponentNames)
            transferFluxStruct(i).(flowComponentNamesNormalized(k)) = transferFlux(i,k);
        end
    end
    transferFlux = transferFluxStruct;

    for k=1:length(flowComponentNames)
        ddtStruct.(flowComponentNamesNormalized(k)) = ddt(k);
    end
    ddt = ddtStruct;

    for k=1:length(flowComponentNames)
        energyStruct.(flowComponentNamesNormalized(k)) = energy(k);
    end
    energy = energyStruct;

    if options.outputFormat == "table"
        forcingFlux = struct2table(forcingFlux);
        forcingFlux.Properties.RowNames = forcingFlux.fancyName;
        forcingFlux.name = [];
        forcingFlux.fancyName = [];

        transferFlux = struct2table(transferFlux);
        transferFlux.Properties.RowNames = transferFlux.fancyName;
        transferFlux.name = [];
        transferFlux.fancyName = [];

        ddt = struct2table(ddt);
        ddt.Properties.RowNames = "ddt";

        energy = struct2table(energy);
        energy.Properties.RowNames = "energy";
    end
end

% forcingNames = self.forcingNames;



function [Tba, Tca, Tcb] = transfersFromFluxes(dA,dB,dC)
if abs(dC) >= abs(dA) && abs(dC) >= abs(dB)
    Tba = 0;
    Tca = dA;
    Tcb = dB;
elseif abs(dB) >= abs(dA) && abs(dB) >= abs(dC)
    Tba = dA;
    Tca = 0;
    Tcb = -dC;
else
    Tba = -dB;
    Tca = -dC;
    Tcb = 0;
end
end

end