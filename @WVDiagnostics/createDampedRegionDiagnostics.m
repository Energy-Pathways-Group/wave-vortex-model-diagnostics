function createDampedRegionDiagnostics(self,options)
arguments
    self WVDiagnostics
    options.stride = 1
    options.timeIndices
    options.shouldUseExplicitAntialiasing logical = false
end

[fpath,fname,~] = fileparts(self.wvpath);
if isempty(fpath)
    fpath = pwd;
end
if options.shouldUseExplicitAntialiasing
    wwgTriadPath = fullfile(fpath,strcat(fname,"-damped-aa.nc"));
    wvt = self.wvt_aa;
else
    wwgTriadPath = fullfile(fpath,strcat(fname,"-damped.nc"));
    wvt = self.wvt;
end


if exist(wwgTriadPath,"file")
    %%%%%%%%%%%%%%%%%%
    % If existing file
    %%%%%%%%%%%%%%%%%%
    ncfile = NetCDFFile(wwgTriadPath);
    t = ncfile.readVariables("t");
    [found, idx] = ismember(t, self.t_wv);
    if ~all(found)
        error('Some entries of t_diag are not exactly in t_wv.');
    end
    if length(t) > 1
        stride = idx(2)-idx(1);
    else
        stride = options.stride;
    end
    timeIndices = (idx(end)+stride):stride:length(self.t_wv);
    if isfield(options,"timeIndices")
        timeIndices = intersect(options.timeIndices,timeIndices);
    end
    outputIndexOffset = length(t);

    ggw_nodamp_var = ncfile.variableWithName("E0_ggw_nodamp");
    wwg_nodamp_var = ncfile.variableWithName("Epm_wwg_nodamp");
else
    %%%%%%%%%%%%%%%%%%
    % If new file
    %%%%%%%%%%%%%%%%%%
    outputIndexOffset = 0;

    if ~isfield(options,"timeIndices")
        timeIndices = 1:options.stride:length(self.t_wv);
    else
        timeIndices = options.timeIndices;
    end

    ncfile = NetCDFFile(wwgTriadPath);
    ncfile.addDimension("kRadial",wvt.kRadial);
    ncfile.addDimension("j",wvt.j);

    varAnnotation = wvt.propertyAnnotationWithName('t');
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    varAnnotation.attributes('standard_name') = 'time';
    varAnnotation.attributes('long_name') = 'time';
    varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
    varAnnotation.attributes('axis') = 'T';
    varAnnotation.attributes('calendar') = 'standard';
    ncfile.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

    dimensionNames = ["j", "kRadial", "t"];
    att = containers.Map(KeyType='char',ValueType='any');
    att("description") = "energy flux of the [g \nabla g]_w triad component for geostrophic modes outside the damping region";
    ggw_nodamp_var = ncfile.addVariable("E0_ggw_nodamp",dimensionNames,type="double",isComplex=false,attributes=att);

    att("description") = "energy flux of the [w \nabla w]_g triad component for wave modes outside the damping region";
    wwg_nodamp_var = ncfile.addVariable("Epm_wwg_nodamp",dimensionNames,type="double",isComplex=false,attributes=att);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the wavenumber/mode mask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

svv = wvt.forcingWithName("adaptive damping");
NoDampMask = (wvt.Kh < (svv.k_damp+svv.k_no_damp)/2) & (wvt.J < (svv.j_damp + svv.j_no_damp)/2);
waveComponent = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
waveComponent.maskAp = waveComponent.maskAp & NoDampMask;
waveComponent.maskAm = waveComponent.maskAm & NoDampMask;

geoComponent = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
geoComponent.maskA0 = geoComponent.maskA0 & NoDampMask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute the WWG triad for %d time indices.\n",length(timeIndices));
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end


    outputIndex = timeIndex + outputIndexOffset;
    self.iTime = timeIndices(timeIndex);
    ncfile.variableWithName('t').setValueAlongDimensionAtIndex(wvt.t,'t',outputIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(waveComponent,waveComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    geo_jk = wvt.transformToRadialWavenumber(E0);
    wwg_nodamp_var.setValueAlongDimensionAtIndex(geo_jk,'t',outputIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(geoComponent,geoComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    wave_jk = wvt.transformToRadialWavenumber(Ep + Em);
    ggw_nodamp_var.setValueAlongDimensionAtIndex(wave_jk,'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end