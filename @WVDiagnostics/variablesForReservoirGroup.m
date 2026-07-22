function [triadVar, forcingVar, energyVar, didAddOutputVariable] = variablesForReservoirGroup(self,options)
% Variables For Reservoir Group.
%
% variablesForReservoirGroup is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Internal — Support functions for createReservoirGroup
% - Declaration: [triadVar, forcingVar, energyVar, didAddOutputVariable] = variablesForReservoirGroup(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter outputfile: file path or name
% - Parameter name: (optional) input argument `name` (default: "reservoir-damped-wave-geo")
% - Parameter flowComponents: input argument `flowComponents`
% - Returns triadVar: output value `triadVar`
% - Returns forcingVar: output value `forcingVar`
% - Returns energyVar: output value `energyVar`
% - Returns didAddOutputVariable: true when the existing group schema required a new output variable
arguments
    self WVDiagnostics
    options.outputfile NetCDFGroup
    options.name string = "reservoir-damped-wave-geo"
    options.flowComponents WVFlowComponent
end

flowComponents = options.flowComponents;
nFlowComponents = length(flowComponents);
triadName = strings(1,nFlowComponents*nFlowComponents*(nFlowComponents+1)/2);
iTriad = 0;
for i=1:nFlowComponents
    for j=1:i
        for k=1:nFlowComponents
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
didAddOutputVariable = false;
if options.outputfile.hasGroupWithName(groupName)
    group = options.outputfile.groupWithName(groupName);
    expectedFlowComponentNames = reshape(string([flowComponents.name]),[],1);
    if ~group.hasAttributeWithName("flow-components")
        error("WVDiagnostics:ReservoirGroupDefinitionMismatch","The existing reservoir group %s has no flow-components definition. Use a different group name or recreate the group.",groupName);
    end
    storedFlowComponentNames = reshape(string(group.attributes("flow-components")),[],1);
    if ~isequal(storedFlowComponentNames,expectedFlowComponentNames)
        error("WVDiagnostics:ReservoirGroupDefinitionMismatch","The existing reservoir group %s defines flow components [%s], but the request defines [%s]. Use a different group name or recreate the group.",groupName,strjoin(storedFlowComponentNames,", "),strjoin(expectedFlowComponentNames,", "));
    end
else
    group = options.outputfile.addGroup(groupName);
    group.addAttribute('flow-components',reshape([flowComponents.name],[],1));
end

for iTriad=1:length(triadName)
    if group.hasVariableWithName(triadName(iTriad))
        triadVar{triadName(iTriad)} = group.variableWithName(triadName(iTriad));
    else
        triadVar{triadName(iTriad)} = group.addVariable(triadName(iTriad),"t",type="double",isComplex=false);
        didAddOutputVariable = true;
    end
end

for i=1:length(forcingNames)
    name = replace(forcingNames(i),"-","_");
    name = replace(name," ","_");
    for k=1:nFlowComponents
        if group.hasVariableWithName(name+"_"+k)
            forcingVar{forcingNames(i)+"_"+k} = group.variableWithName(name+"_"+k);
        else
            forcingVar{forcingNames(i)+"_"+k} = group.addVariable(name+"_"+k,"t",type="double",isComplex=false);
            didAddOutputVariable = true;
        end
    end
end

for k=1:nFlowComponents
    if group.hasVariableWithName("E_"+k)
        energyVar{"E_"+k} = group.variableWithName("E_"+k);
    else
        energyVar{"E_"+k} = group.addVariable("E_"+k,"t",type="double",isComplex=false);
        didAddOutputVariable = true;
    end
end
end
