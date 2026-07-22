function createGeostrophicFluxGroup(self,options)
% Create a geostrophic flux group in the diagnostics NetCDF and populate it.
%
% Computes the geostrophic fluxes from geostrophicFlux, valid only for
% constant stratification. Existing completed time slices are preserved by
% default, and incomplete slices are computed or repaired.
%
% - Topic: Diagnostics Generation
% - Declaration: createGeostrophicFluxGroup(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter outputfile: NetCDFGroup to write into (default: self.diagfile).
% - Parameter timeIndices: Time indices to process (default: all times in diagnostics file).
% - Parameter shouldOverwriteExisting: Recompute requested completed slices when true (default: false).
arguments
    self WVDiagnostics
    options.outputfile NetCDFGroup
    options.timeIndices
    options.shouldOverwriteExisting (1,1) logical = false
end

if self.diagnosticsHasExplicitAntialiasing
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end

if ~isfield(options,"outputfile")
    options.outputfile = self.diagfile;
end

isNewGroup = ~options.outputfile.hasGroupWithName("geostrophic-flux");
if isNewGroup
    group = options.outputfile.addGroup("geostrophic-flux");
else
    group = options.outputfile.groupWithName("geostrophic-flux");
end

dimensionNames = ["j", "kRadial", "t"];
variableNames = ["ggg", "ggw", "ggw_tx", "wwg_tx"];
outputVariables = cell(size(variableNames));
didAddOutputVariable = false;
for iVariable = 1:length(variableNames)
    if group.hasVariableWithName(variableNames(iVariable))
        outputVariables{iVariable} = group.variableWithName(variableNames(iVariable));
    else
        outputVariables{iVariable} = group.addVariable(variableNames(iVariable),dimensionNames,type="double",isComplex=false);
        didAddOutputVariable = true;
    end
end
[ggg,ggw,ggw_tx,wwg_tx] = outputVariables{:};

diagnosticTimes = self.diagfile.readVariables("t");
if ~isfield(options,"timeIndices")
    requestedTimeIndices = 1:length(diagnosticTimes);
else
    requestedTimeIndices = options.timeIndices;
end
[timeIndices,completionVariable] = self.timeIndicesToComputeForDiagnosticGroup(group,outputVariables,requestedTimeIndices,shouldOverwriteExisting=options.shouldOverwriteExisting,isNewGroup=isNewGroup,didAddOutputVariable=didAddOutputVariable);
if isempty(timeIndices)
    fprintf("Geostrophic flux diagnostics are already complete for the requested time indices.\n");
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute geostrophic fluxes for %d time indices.\n",length(timeIndices));
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex + 1);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndices(timeIndex);
    self.iTime = timeIndices(timeIndex);

    [~, spectralFlux] = self.geostrophicFlux();
    ggg.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggg),'t',outputIndex);
    ggw.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggw),'t',outputIndex);
    ggw_tx.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggw_tx),'t',outputIndex);
    wwg_tx.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.wwg_tx),'t',outputIndex);
    completionVariable.setValueAlongDimensionAtIndex(diagnosticTimes(outputIndex),'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end
