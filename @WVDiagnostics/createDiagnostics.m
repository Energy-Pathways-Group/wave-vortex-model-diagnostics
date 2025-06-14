function createDiagnostics(self,options)
arguments
    self WVDiagnostics 
    options.stride = 1
    options.timeIndices
    options.outpath
    options.shouldMeasureAntialiasingFlux logical = false
    options.shouldUseHigherOrderFlux logical = false
end

if exist(self.diagfile,"file")
    warning('A diagnostics file already exists. Returning.');
    return
end

if options.shouldMeasureAntialiasingFlux
    [wvt_lowres, ncfile] = WVTransform.waveVortexTransformFromFile(self.wvfile.path,iTime=Inf);
    wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
else
    [wvt, ncfile] = WVTransform.waveVortexTransformFromFile(self.wvfile.path,iTime=Inf);
end
tDim = ncfile.readVariables('wave-vortex/t');
if ~isfield(options,"timeIndices")
    timeIndices = 1:options.stride:length(tDim);
else
    timeIndices = options.timeIndices;
end

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APEOperation(wvt));
wvt.addOperation(APVOperation());
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

%% setup diagnostic output file
diagfile = NetCDFFile(self.diagpath);
self.diagfile = diagfile;

dimensionNames = ["j","kRadial"];
for iDim=1:length(dimensionNames)
    dimAnnotation = wvt.dimensionAnnotationWithName(dimensionNames(iDim));
    dimAnnotation.attributes('units') = dimAnnotation.units;
    dimAnnotation.attributes('long_name') = dimAnnotation.description;
    diagfile.addDimension(dimAnnotation.name,wvt.(dimAnnotation.name),attributes=dimAnnotation.attributes);
end

varAnnotation = wvt.propertyAnnotationWithName('t');
varAnnotation.attributes('units') = varAnnotation.units;
varAnnotation.attributes('long_name') = varAnnotation.description;
varAnnotation.attributes('standard_name') = 'time';
varAnnotation.attributes('long_name') = 'time';
varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
varAnnotation.attributes('axis') = 'T';
varAnnotation.attributes('calendar') = 'standard';
diagfile.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

% 1. Measures of energy, APV and enstrophy
EnergyByComponent = configureDictionary("string","cell");
flowComponentNames = wvt.flowComponentNames;
for i=1:length(flowComponentNames)
    comp = wvt.flowComponentWithName(flowComponentNames(i));
    name = "E_" + comp.abbreviatedName;
    EnergyByComponent{name} = diagfile.addVariable(name,"t",type="double",isComplex=false);
end
EnergyByComponent{"KE_g"} = diagfile.addVariable("KE_g","t",type="double",isComplex=false);
EnergyByComponent{"PE_g"} = diagfile.addVariable("PE_g","t",type="double",isComplex=false);
variable_ke = diagfile.addVariable("ke","t",type="double",isComplex=false);
variable_pe = diagfile.addVariable("pe_quadratic","t",type="double",isComplex=false);
variable_ape = diagfile.addVariable("ape","t",type="double",isComplex=false);
variable_z = diagfile.addVariable("enstrophy_quadratic","t",type="double",isComplex=false);
variable_apv = diagfile.addVariable("enstrophy_apv","t",type="double",isComplex=false);

% 2. Forcing fluxes
EnergyFlux = configureDictionary("string","cell");
EnstrophyFlux = configureDictionary("string","cell");
EnergyFluxTrue = configureDictionary("string","cell");
forcingNames = wvt.forcingNames;
dimensionNames = ["j", "kRadial", "t"];
for i=1:length(forcingNames)
    name = replace(forcingNames(i),"-","_");
    name = replace(name," ","_");
    EnergyFlux{forcingNames(i)}.Ep = diagfile.addVariable("Ep_" + name,dimensionNames,type="double",isComplex=false);
    EnergyFlux{forcingNames(i)}.Em = diagfile.addVariable("Em_" + name,dimensionNames,type="double",isComplex=false);
    EnergyFlux{forcingNames(i)}.KE0 = diagfile.addVariable("KE0_" + name,dimensionNames,type="double",isComplex=false);
    EnergyFlux{forcingNames(i)}.PE0 = diagfile.addVariable("PE0_" + name,dimensionNames,type="double",isComplex=false);
    EnstrophyFlux{forcingNames(i)} = diagfile.addVariable("Z0_" + name,dimensionNames,type="double",isComplex=false);
    EnergyFluxTrue{forcingNames(i)} = diagfile.addVariable("E_" + name,"t",type="double",isComplex=false);
end

% 3. Triads
triadFlowComponents = [wvt.flowComponentWithName('wave'); wvt.flowComponentWithName('inertial'); wvt.flowComponentWithName('geostrophic'); wvt.flowComponentWithName('mda')];
EnergyTriads = configureDictionary("string","cell");
EnstrophyTriads = configureDictionary("string","cell");
dimensionNames = ["j", "kRadial", "t"];
for i=1:length(triadFlowComponents)
    for j=1:length(triadFlowComponents)
        name = triadFlowComponents(i).abbreviatedName + "_" + triadFlowComponents(j).abbreviatedName;
        EnergyTriads{name}.Ep = diagfile.addVariable("Ep_" + name,dimensionNames,type="double",isComplex=false);
        EnergyTriads{name}.Em = diagfile.addVariable("Em_" + name,dimensionNames,type="double",isComplex=false);
        EnergyTriads{name}.KE0 = diagfile.addVariable("KE0_" + name,dimensionNames,type="double",isComplex=false);
        EnergyTriads{name}.PE0 = diagfile.addVariable("PE0_" + name,dimensionNames,type="double",isComplex=false);
        EnstrophyTriads{name} = diagfile.addVariable("Z0_" + name,dimensionNames,type="double",isComplex=false);
    end
end

%% loop through time computing diagnostics
for outputIndex = 1:length(timeIndices)
    if mod(outputIndex,10) == 1
        fprintf("%d..",outputIndex);
    end
    iTime = timeIndices(outputIndex);
    if options.shouldMeasureAntialiasingFlux
        wvt_lowres.initFromNetCDFFile(ncfile,iTime=iTime);
        [wvt.A0,wvt.Ap,wvt.Am] = wvt_lowres.spectralVariableWithResolution(wvt,wvt_lowres.A0,wvt_lowres.Ap,wvt_lowres.Am);
    else
        wvt.initFromNetCDFFile(ncfile,iTime=iTime);
    end

    diagfile.variableWithName('t').setValueAlongDimensionAtIndex(wvt.t,'t',outputIndex);

    % 1. Measures of energy, APV and enstrophy
    for i=1:length(flowComponentNames)
        comp = wvt.flowComponentWithName(flowComponentNames(i));
        name = "E_" + comp.abbreviatedName;
        E = wvt.totalEnergyOfFlowComponent(comp);
        EnergyByComponent{name}.setValueAlongDimensionAtIndex(E,'t',outputIndex);
    end
    EnergyByComponent{"KE_g"}.setValueAlongDimensionAtIndex(wvt.geostrophicKineticEnergy,'t',outputIndex);
    EnergyByComponent{"PE_g"}.setValueAlongDimensionAtIndex(wvt.geostrophicPotentialEnergy,'t',outputIndex);
    [u,v,ape,apv] = wvt.variableWithName('u','v','ape','apv');
    ke = (u.^2 + v.^2)/2;
    variable_ke.setValueAlongDimensionAtIndex(int_vol(ke),'t',outputIndex);
    variable_pe.setValueAlongDimensionAtIndex(int_vol(shiftdim(wvt.N2,-2).*wvt.eta.*wvt.eta/2),'t',outputIndex);
    variable_ape.setValueAlongDimensionAtIndex(int_vol(ape),'t',outputIndex);
    variable_z.setValueAlongDimensionAtIndex(wvt.totalEnstrophy,'t',outputIndex);
    variable_apv.setValueAlongDimensionAtIndex(0.5*int_vol(apv.*apv),'t',outputIndex);

    % 2. Forcing fluxes
    F = wvt.fluxForForcing();
    for i=1:length(forcingNames)
        [Ep,Em,KE0,PE0] = wvt.energyFluxFromNonlinearFlux(F{forcingNames(i)}.Fp,F{forcingNames(i)}.Fm,F{forcingNames(i)}.F0);
        [Ep_jk,Em_jk,KE0_jk,PE0_jk] = wvt.transformToRadialWavenumber(Ep,Em,KE0,PE0);
        if isa(wvt,"WVTransformHydrostatic")
            [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(i));
            F_density = wvt.u .* Fu + wvt.v .* Fv+ wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;
        elseif isa(wvt,"WVTransformBoussinesq")
            [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(i));
            F_density = wvt.u .* Fu + wvt.v .* Fv +  wvt.w .* Fw + wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;
        else
            error("Transform not yet supported.");
        end

        if forcingNames(i) == "nonlinear advection"
            F_density = F_density + wvt.w .* shiftdim(wvt.N2,-2) .* (wvt.eta_true-wvt.eta);
        end

        EnergyFluxTrue{forcingNames(i)}.setValueAlongDimensionAtIndex(int_vol(F_density),'t',outputIndex);

        EnergyFlux{forcingNames(i)}.Ep.setValueAlongDimensionAtIndex(Ep_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.Em.setValueAlongDimensionAtIndex(Em_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.KE0.setValueAlongDimensionAtIndex(KE0_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.PE0.setValueAlongDimensionAtIndex(PE0_jk,'t',outputIndex);

        Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(i)}.F0 .* conj(wvt.A0) );
        EnstrophyFlux{forcingNames(i)}.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(Z),'t',outputIndex);
    end

    % 3. Triads
    for i=1:length(triadFlowComponents)
        for j=1:length(triadFlowComponents)
            if options.shouldUseHigherOrderFlux
                [Fp,Fm,F0] = wvt.rk4NonlinearFluxForFlowComponents(triadFlowComponents(i),triadFlowComponents(j));
            else
                [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(triadFlowComponents(i),triadFlowComponents(j));
            end
            [Ep,Em,KE0,PE0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
            Z0 = wvt.enstrophyFluxFromNonlinearFlux(F0);
            [Ep_jk,Em_jk,KE0_jk,PE0_jk,Z0_jk] = wvt.transformToRadialWavenumber(Ep,Em,KE0,PE0,Z0);

            name = triadFlowComponents(i).abbreviatedName + "_" + triadFlowComponents(j).abbreviatedName;
            EnergyTriads{name}.Ep.setValueAlongDimensionAtIndex(Ep_jk,'t',outputIndex);
            EnergyTriads{name}.Em.setValueAlongDimensionAtIndex(Em_jk,'t',outputIndex);
            EnergyTriads{name}.KE0.setValueAlongDimensionAtIndex(KE0_jk,'t',outputIndex);
            EnergyTriads{name}.PE0.setValueAlongDimensionAtIndex(PE0_jk,'t',outputIndex);
            EnstrophyTriads{name}.setValueAlongDimensionAtIndex(Z0_jk,'t',outputIndex);
        end
    end
end
fprintf("\n");
diagfile.close();
end