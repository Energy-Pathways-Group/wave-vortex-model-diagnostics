basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% wvd1 = WVDiagnostics(basedir + replace(getRunParameters(1),"256","512") + ".nc");
wvd9 = WVDiagnostics(basedir + replace(getRunParameters(9),"256","256") + ".nc");
% wvd18 = WVDiagnostics(basedir + replace(getRunParameters(18),"256","256") + ".nc");

%%
wvd = wvd9;
wvt = wvd.wvt;

timeIndices = 2751:3001;

diagfile = wvd.diagfile;
if ~diagfile.hasVariableWithName("T_ggw_gw")
    T_ggw_gw = diagfile.addVariable("T_ggw_gw","t",type="double",isComplex=false);
    T_ggg_gg = diagfile.addVariable("T_ggg_gg","t",type="double",isComplex=false);
    T_ggd_gd = diagfile.addVariable("T_ggd_gd","t",type="double",isComplex=false);
    T_www_ww = diagfile.addVariable("T_www_ww","t",type="double",isComplex=false);
    T_wwg_wg = diagfile.addVariable("T_wwg_wg","t",type="double",isComplex=false);
    T_wwd_wd = diagfile.addVariable("T_wwd_wd","t",type="double",isComplex=false);
    T_ddw_dw = diagfile.addVariable("T_ddw_dw","t",type="double",isComplex=false);
    T_ddg_dg = diagfile.addVariable("T_ddg_dg","t",type="double",isComplex=false);
    T_ddd_dd = diagfile.addVariable("T_ddd_dd","t",type="double",isComplex=false);
    T_gwd_gw = diagfile.addVariable("T_gwd_gw","t",type="double",isComplex=false);
    T_gwd_gd = diagfile.addVariable("T_gwd_gd","t",type="double",isComplex=false);
    T_gwd_wd = diagfile.addVariable("T_gwd_wd","t",type="double",isComplex=false);
else
    T_ggw_gw = diagfile.variableWithName("T_ggw_gw");
    T_ggg_gg = diagfile.variableWithName("T_ggg_gg");
    T_ggd_gd = diagfile.variableWithName("T_ggd_gd");
    T_www_ww = diagfile.variableWithName("T_www_ww");
    T_wwg_wg = diagfile.variableWithName("T_wwg_wg");
    T_wwd_wd = diagfile.variableWithName("T_wwd_wd");
    T_ddw_dw = diagfile.variableWithName("T_ddw_dw");
    T_ddg_dg = diagfile.variableWithName("T_ddg_dg");
    T_ddd_dd = diagfile.variableWithName("T_ddd_dd");
    T_gwd_gw = diagfile.variableWithName("T_gwd_gw");
    T_gwd_gd = diagfile.variableWithName("T_gwd_gd");
    T_gwd_wd = diagfile.variableWithName("T_gwd_wd");
end

%%

svv = wvt.forcingWithName("adaptive damping");
NoDampMask = (wvt.Kh < (svv.k_damp+svv.k_no_damp)/2) & (wvt.J < (svv.j_damp + svv.j_no_damp)/2);

waveComponent = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
waveComponent.maskAp = waveComponent.maskAp & NoDampMask;
waveComponent.maskAm = waveComponent.maskAm & NoDampMask;

geoComponent = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
geoComponent.maskA0 = geoComponent.maskA0 & NoDampMask;

dampComponent = wvt.totalFlowComponent;
dampComponent.maskAp = dampComponent.maskAp & ~NoDampMask;
dampComponent.maskAm = dampComponent.maskAm & ~NoDampMask;
dampComponent.maskA0 = dampComponent.maskA0 & ~NoDampMask;

waveEnergy = @(Ep,Em,E0) sum(waveComponent.maskAp(:).*Ep(:) + waveComponent.maskAm(:).*Em(:));
geoEnergy = @(Ep,Em,E0) sum(geoComponent.maskA0(:).*E0(:));
dampEnergy = @(Ep,Em,E0) sum(dampComponent.maskAp(:).*Ep(:) + dampComponent.maskAm(:).*Em(:) + dampComponent.maskA0(:).*E0(:));

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
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex + 1);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndices(timeIndex);
    wvd.iTime = timeIndices(timeIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(geoComponent,geoComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    T_ggw_gw.setValueAlongDimensionAtIndex(waveEnergy(Ep,Em,E0),'t',outputIndex);
    T_ggg_gg.setValueAlongDimensionAtIndex(geoEnergy(Ep,Em,E0),'t',outputIndex);
    T_ggd_gd.setValueAlongDimensionAtIndex(dampEnergy(Ep,Em,E0),'t',outputIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(waveComponent,waveComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    T_www_ww.setValueAlongDimensionAtIndex(waveEnergy(Ep,Em,E0),'t',outputIndex);
    T_wwg_wg.setValueAlongDimensionAtIndex(geoEnergy(Ep,Em,E0),'t',outputIndex);
    T_wwd_wd.setValueAlongDimensionAtIndex(dampEnergy(Ep,Em,E0),'t',outputIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(dampComponent,dampComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    T_ddw_dw.setValueAlongDimensionAtIndex(waveEnergy(Ep,Em,E0),'t',outputIndex);
    T_ddg_dg.setValueAlongDimensionAtIndex(geoEnergy(Ep,Em,E0),'t',outputIndex);
    T_ddd_dd.setValueAlongDimensionAtIndex(dampEnergy(Ep,Em,E0),'t',outputIndex);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(dampComponent,waveComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEg_dwg = geoEnergy(Ep,Em,E0);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(waveComponent,dampComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEg_wdg = geoEnergy(Ep,Em,E0);
    dEg = dEg_dwg + dEg_wdg;

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(dampComponent,geoComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEw_dgw = waveEnergy(Ep,Em,E0);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(geoComponent,dampComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEw_gdw = waveEnergy(Ep,Em,E0);
    dEw = dEw_dgw + dEw_gdw;

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(geoComponent,waveComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEd_gwd = dampEnergy(Ep,Em,E0);

    [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(waveComponent,geoComponent);
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
    dEd_wgd = dampEnergy(Ep,Em,E0);
    dEd = dEd_gwd + dEd_wgd;

    T_gwd_gw.setValueAlongDimensionAtIndex((dEg - dEw)/3,'t',outputIndex);
    T_gwd_gd.setValueAlongDimensionAtIndex((dEg - dEd)/3,'t',outputIndex);
    T_gwd_wd.setValueAlongDimensionAtIndex((dEw - dEd)/3,'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

%%
T_wg = T_wwg_wg.value - T_ggw_gw.value + T_gwd_gw.value;
T_dg = T_ddg_dg.value - T_ggd_gd.value + T_gwd_gd.value;
T_dw = T_ddw_dw.value - T_wwd_wd.value + T_gwd_wd.value;

mean(T_wg(timeIndices))/wvd.flux_scale
mean(T_dg(timeIndices))/wvd.flux_scale
mean(T_dw(timeIndices))/wvd.flux_scale