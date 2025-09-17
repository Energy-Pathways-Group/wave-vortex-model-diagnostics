function createWWGTriadDiagnostic(self,options)
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
    wwgTriadPath = fullfile(fpath,strcat(fname,"-wwg-aa.nc"));
    wvt = self.wvt_aa;
else
    wwgTriadPath = fullfile(fpath,strcat(fname,"-wwg.nc"));
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
        stride = didx(1);
    else
        stride = options.stride;
    end
    timeIndices = (idx(end)+stride):stride:length(self.t_wv);
    if isfield(options,"timeIndices")
        timeIndices = intersect(options.timeIndices,timeIndices);
    end
    outputIndexOffset = length(t);

    ww_var = ncfile.variableWithName("wavewave_jk");
    geo_var = ncfile.variableWithName("geo_jk");
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
    ww_var = ncfile.addVariable("wavewave_jk",dimensionNames,type="double",isComplex=false);
    geo_var = ncfile.addVariable("geo_jk",dimensionNames,type="double",isComplex=false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the wavenumber/mode mask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = wvt.kRadial;
dk = k(2)-k(1);

nK = length(k);
Kh = wvt.Kh;

maskK = false(nK,wvt.Nkl);
for iK = 1:1:nK
    indicesForK = Kh(1,:) < k(iK)+dk/2;
    maskK(iK,indicesForK) = true;
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
fprintf("Starting loop to compute the WWG triad for %d time indices, over Nj=%d, Nk=%d.\n",length(timeIndices),wvt.Nj,nK);
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

    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1,1);
    geo_jk = wvt.transformToRadialWavenumber(E0);
    geo_var.setValueAlongDimensionAtIndex(geo_jk,'t',outputIndex);

    wavewave = zeros(wvt.Nj,nK);

    fprintf("%d: iTime=%d. k=",timeIndex,self.iTime);
    for i=1:nK
        if mod(i,10) == 0
            fprintf("%d..",i)
        end
        for j=1:wvt.Nj
            E0 = WVDiagnostics.waveWaveGeostrophicEnergyForMode(wvt,maskK(i,:),maskK(i,:),j);
            wavewave(j,i) = sum(E0(:));
        end
    end
    fprintf("\n");
    
    ww_var.setValueAlongDimensionAtIndex(wavewave,'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end