function [triadVar, forcingVar, energyVar] = variablesForReservoirGroup(self,options)
arguments
    self WVDiagnostics
    options.outputfile NetCDFGroup
    options.name string = "reservoir-damped-wave-geo"
    options.flowComponents WVFlowComponent
end

flowComponents = options.flowComponents;
iTriad = 0;
for i=1:length(flowComponents)
    for j=1:i
        for k=1:length(flowComponents)
            iTriad = iTriad + 1;
            triadName(iTriad) = "T_" + i + "_" + j + "_" + k;
        end
    end
end

if self.diagnosticsHasExplicitAntialiasing
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end

groupName = options.name;
forcingNames = wvt.forcingNames;
triadVar = configureDictionary("string","cell");
forcingVar = configureDictionary("string","cell");
energyVar = configureDictionary("string","cell");
if options.outputfile.hasGroupWithName(groupName)
    group = options.outputfile.groupWithName(groupName);
    for iTriad=1:length(triadName)
        triadVar{triadName(iTriad)} = group.variableWithName(triadName(iTriad));
    end
    for i=1:length(forcingNames)
        name = replace(forcingNames(i),"-","_");
        name = replace(name," ","_");
        for k=1:length(flowComponents)
            forcingVar{forcingNames(i)+"_"+k} = group.variableWithName(name+"_"+k);
        end
    end
    for k=1:length(flowComponents)
        energyVar{"E_"+k} = group.variableWithName("E_"+k);
    end
else
    group = options.outputfile.addGroup(groupName);
    group.addAttribute('flow-components',reshape([flowComponents.name],[],1));

    % Add a variable for each triad component
    for iTriad=1:length(triadName)
        triadVar{triadName(iTriad)} = group.addVariable(triadName(iTriad),"t",type="double",isComplex=false);
    end

    % Add a variable for each forcing for each reservoir
    for i=1:length(forcingNames)
        name = replace(forcingNames(i),"-","_");
        name = replace(name," ","_");
        for k=1:length(flowComponents)
            forcingVar{forcingNames(i)+"_"+k} = group.addVariable(name+"_"+k,"t",type="double",isComplex=false);
        end
    end

    for k=1:length(flowComponents)
        energyVar{"E_"+k} = group.addVariable("E_"+k,"t",type="double",isComplex=false);
    end
end
end