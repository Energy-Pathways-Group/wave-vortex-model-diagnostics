function create1DMirrorFluxes(self,options)
arguments
    self WVDiagnostics
    options.stride = 1
    options.timeIndices
end

if ~exist(self.diagpath,"file")
    error("No existing diagnostics file found.");
end
diagfile = self.diagfile;
wvt = self.wvt;

[kp,bins_0,bins_pm] = self.sparsePseudoRadialAxis;

if diagfile.hasGroupWithName("mirror-fluxes-1d")
    %%%%%%%%%%%%%%%%%%
    % If existing file
    %%%%%%%%%%%%%%%%%%
    group = diagfile.groupWithName("mirror-fluxes-1d");
    t = group.readVariables("t");
    [found, idx] = ismember(t, self.t_wv);
    if ~all(found)
        error('Some entries of mirror-fluxes-1d/t are not exactly in t_wv.');
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

    pi_w_wwg_kp = group.variableWithName("pi_w_wwg_kp");
    F_wwg_kp = group.variableWithName("F_wwg_kp");
    pi_g_ggw_kp = group.variableWithName("pi_g_ggw_kp");
    F_ggw_kp = group.variableWithName("F_ggw_kp");
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

    group = diagfile.addGroup("mirror-fluxes-1d");
    group.addDimension("kp",kp);

    varAnnotation = wvt.propertyAnnotationWithName('t');
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    varAnnotation.attributes('standard_name') = 'time';
    varAnnotation.attributes('long_name') = 'time';
    varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
    varAnnotation.attributes('axis') = 'T';
    varAnnotation.attributes('calendar') = 'standard';
    group.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

    dimensionNames = ["kp", "t"];
    pi_w_wwg_kp = group.addVariable("pi_w_wwg_kp",dimensionNames,type="double",isComplex=false);
    F_wwg_kp = group.addVariable("F_wwg_kp",dimensionNames,type="double",isComplex=false);
    pi_g_ggw_kp = group.addVariable("pi_g_ggw_kp",dimensionNames,type="double",isComplex=false);
    F_ggw_kp = group.addVariable("F_ggw_kp",dimensionNames,type="double",isComplex=false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the wavenumber/mode mask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valid = ~isnan(bins_0);
S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

valid = ~isnan(bins_pm);
S_pm = sparse(find(valid), bins_pm(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

mask_0 = false(wvd.wvt.Nj,wvt.Nkl,length(kp));
mask_pm = false(wvd.wvt.Nj,wvt.Nkl,length(kp));
for iK = 1:1:length(kp)
    mask_0(:,:,iK) = (bins_0 <= iK);
    mask_pm(:,:,iK) = (bins_pm <= iK);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute the 1d mirror fluxes for %d time indices, over N=%d, Nk=%d.\n",length(timeIndices),length(kp));
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
    group.variableWithName('t').setValueAlongDimensionAtIndex(wvt.t,'t',outputIndex);

    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1);
    F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
    F_wwg_kp.setValueAlongDimensionAtIndex(F_wwg_kp_val,'t',outputIndex);

    Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,1);
    F_ggw_kp_val = reshape(Epm(:).' * S_pm,[],1);
    F_ggw_kp.setValueAlongDimensionAtIndex(F_ggw_kp_val,'t',outputIndex);

    pi_w_wwg_kp_val = zeros(length(kp),1);
    pi_g_ggw_kp_val = zeros(length(kp),1);

    for i=1:length(kp)
        E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_pm(:,:,i));
        pi_w_wwg_kp_val(i) = sum(E0(:));
        Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_0(:,:,i));
        pi_g_ggw_kp_val(i) = sum(Epm(:));
    end
    pi_w_wwg_kp.setValueAlongDimensionAtIndex(pi_w_wwg_kp_val,'t',outputIndex);
    pi_g_ggw_kp.setValueAlongDimensionAtIndex(pi_g_ggw_kp_val,'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end