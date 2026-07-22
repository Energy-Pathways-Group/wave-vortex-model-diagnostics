function [timeIndices, completionVariable] = timeIndicesToComputeForDiagnosticGroup(self,group,outputVariables,requestedTimeIndices,options)
% Select incomplete diagnostic time indices and maintain completion metadata.
arguments
    self WVDiagnostics
    group NetCDFGroup
    outputVariables cell
    requestedTimeIndices double
    options.shouldOverwriteExisting (1,1) logical = false
    options.isNewGroup (1,1) logical = false
    options.didAddOutputVariable (1,1) logical = false
end

diagnosticTimes = self.diagfile.readVariables("t");
diagnosticTimes = diagnosticTimes(:);
requestedTimeIndices = reshape(requestedTimeIndices,1,[]);

if ~isempty(requestedTimeIndices)
    isValid = isfinite(requestedTimeIndices) & requestedTimeIndices == fix(requestedTimeIndices) & requestedTimeIndices >= 1 & requestedTimeIndices <= length(diagnosticTimes);
    if ~all(isValid)
        error("WVDiagnostics:InvalidTimeIndices","Time indices must be finite positive integers between 1 and %d.",length(diagnosticTimes));
    end
end
requestedTimeIndices = unique(requestedTimeIndices,"stable");

hasCompletionVariable = group.hasVariableWithName("t");
if hasCompletionVariable
    completionVariable = group.variableWithName("t");
else
    timeVariable = self.diagfile.variableWithName("t");
    completionVariable = group.addVariable("t","t",type="double",isComplex=false,attributes=timeVariable.attributes);
end

fillValue = netcdf.getConstant('NC_FILL_DOUBLE');
if hasCompletionVariable && options.didAddOutputVariable
    completionVariable.setValueAlongDimensionAtIndex(fillValue*ones(size(diagnosticTimes)),"t",1:length(diagnosticTimes));
elseif ~hasCompletionVariable && ~options.isNewGroup && ~options.didAddOutputVariable
    for outputIndex = 1:length(diagnosticTimes)
        isComplete = true;
        for iVariable = 1:length(outputVariables)
            value = outputVariables{iVariable}.valueAlongDimensionAtIndex("t",outputIndex);
            if isempty(value) || any(value(:) == fillValue)
                isComplete = false;
                break
            end
        end
        if isComplete
            completionVariable.setValueAlongDimensionAtIndex(diagnosticTimes(outputIndex),"t",outputIndex);
        end
    end
end

if options.shouldOverwriteExisting
    timeIndices = requestedTimeIndices;
    return
end

completionTimes = completionVariable.value;
completionTimes = completionTimes(:);
completed = false(size(diagnosticTimes));
nCompletionTimes = min(length(completionTimes),length(diagnosticTimes));
completed(1:nCompletionTimes) = completionTimes(1:nCompletionTimes) == diagnosticTimes(1:nCompletionTimes);
timeIndices = requestedTimeIndices(~completed(requestedTimeIndices));
end
