classdef WVDiagnostics < handle
    %WVDiagnostics Produces diagnostics and figures from WVModel output
    %   This is a collection of diagnostic tools for analyzing model output
    properties
        wvfile
        diagfile
        wvt
    end

    properties (Dependent)
        t_diag
        t_wv
    end

    properties (SetObservable)
        iTime = Inf
    end

    methods
        function self = WVDiagnostics(filename)
            arguments
                filename
            end
            [self.wvt, self.wvfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=Inf);
            [fpath,fname,~] = fileparts(filename);
            if ~isempty(fpath)
                diagfilename = fullfile(fpath,strcat(fname,"-diagnostics.nc"));
            else
                diagfilename = fullfile(strcat(fname,"-diagnostics.nc"));
            end
            if exist(diagfilename,"file")
                self.diagfile = NetCDFFile(diagfilename);
            else
                warning("No diagnostics file found. Some functionality will not be available.")
            end

            addlistener(self,'iTime','PostSet',@WVDiagnostics.iTimeChanged);
        end

        function t = get.t_diag(self)
            t = self.diagfile.readVariables('t');
        end

        function t = get.t_wv(self)
            t = self.wvfile.readVariables('t');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (over time)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fig = plotMooringRotarySpectrum(self)
        fig = plotFluidStateMultipanel(self,options)
        
        function fig = plotEnstrophyOverTime(self,options)
            arguments
                self WVDiagnostics
                options.tscale = 86400;
                options.tscale_units = "days";
                options.visible = "on"
            end
            [Z_quadratic,Z_apv] = self.diagfile.readVariables('enstrophy_quadratic','enstrophy_apv');
            t = self.diagfile.readVariables('t');
            zscale = self.wvt.f^2;
            zscale_units = "f^2";
            fig = figure(Visible=options.visible);
            plot(t/options.tscale,Z_quadratic/zscale,LineWidth=2), hold on
            plot(t/options.tscale,Z_apv/zscale,LineWidth=2)
            legend('quadratic','apv')

            xlabel("time (" + options.tscale_units + ")")
            ylabel("enstrophy (" + zscale_units + ")")
            xlim([min(t) max(t)]/options.tscale);
        end

        function fig = plotEnergyForReservoirOverTime(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te", "te_quadratic"]
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74
                options.escale_units = "GM"
                options.visible = "on"
            end
            reservoirs = self.energyForReservoirOverTime(reservoirNames=options.reservoirNames);
            t = self.diagfile.readVariables('t');

            fig = figure(Visible=options.visible);

            for iReservoir = 1:length(reservoirs)
                switch reservoirs(iReservoir).name
                    case "te"
                        plot(t/options.tscale,reservoirs(iReservoir).energy/options.escale,LineWidth=2, Color=[0 0 0]), hold on
                    case "te_quadratic"
                        plot(t/options.tscale,reservoirs(iReservoir).energy/options.escale,LineWidth=2, Color=[0 0 0], LineStyle="-."), hold on
                    otherwise
                        plot(t/options.tscale,reservoirs(iReservoir).energy/options.escale,LineWidth=2), hold on
                end
            end
            legend(reservoirs.fancyName)

            xlabel("time (" + options.tscale_units + ")")
            ylabel("energy (" + options.escale_units + ")")
            xlim([min(t) max(t)]/options.tscale);
        end

        function fig = plotForcingFluxForReservoirOverTime(self,options)
            % wvd.plotEnergyFluxForReservoirOverTime(filter=@(v) movmean(v,51));
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"]
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74/(86400*365)
                options.escale_units = "GM/yr"
                options.visible = "on"
                options.filter = @(v) v;
            end
            forcing_fluxes = self.forcingFluxesForReservoirOverTime(reservoirNames=options.reservoirNames);
            t = self.diagfile.readVariables('t');

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.reservoirNames),1,TileSpacing="compact");
            for iReservoir = 1:length(options.reservoirNames)
                nexttile(tl);
                for iForce = 1:length(forcing_fluxes)
                    plot(t/options.tscale,options.filter(forcing_fluxes(iForce).(options.reservoirNames(iReservoir))/options.escale)), hold on
                end
                legend(forcing_fluxes.fancyName)

                fancyName = self.fancyNameForName(options.reservoirNames(iReservoir));
                xlabel("time (" + options.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + options.escale_units + ")")
                xlim([min(t) max(t)]/options.tscale);
            end
        end

        function fig = plotInertialFluxForReservoirOverTime(self,options)
            % wvd.plotInertialFluxForReservoirOverTime(filter=@(v) movmean(v,51));
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"]
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74/(86400*365)
                options.escale_units = "GM/yr"
                options.visible = "on"
                options.filter = @(v) v;
            end
            inertial_fluxes = self.inertialFluxesForReservoirOverTime(reservoirNames=options.reservoirNames);
            t = self.diagfile.readVariables('t');

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.reservoirNames),1,TileSpacing="compact");
            for iReservoir = 1:length(options.reservoirNames)
                nexttile(tl);
                for iForce = 1:length(inertial_fluxes)
                    plot(t/options.tscale,options.filter(inertial_fluxes(iForce).(options.reservoirNames(iReservoir))/options.escale)), hold on
                end
                legend(inertial_fluxes.fancyName)

                fancyName = self.fancyNameForName(options.reservoirNames(iReservoir));
                xlabel("time (" + options.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + options.escale_units + ")")
                xlim([min(t) max(t)]/options.tscale);
            end
        end

        function fig = plotSourcesSinksReservoirsDiagram(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"]
                options.flux_scale = 3.74/(86400*365)
                options.flux_scale_units = "GM/yr"
                options.energy_scale = 3.74
                options.energy_scale_units = "GM"
                options.customNames = configureDictionary("string","string")
                options.fluxTolerance = 1e-2;
                options.shouldShowUnits = true;
                options.timeIndices = Inf;
                options.shouldShowReservoirEnergy = true
                options.title = "Energy Pathways";
                options.visible = "on"
            end
            forcing_fluxes = self.forcingFluxesForReservoirSpatialTemporalAverage(reservoirNames=options.reservoirNames,timeIndices=options.timeIndices);

            col = configureDictionary("string","cell");
            col{"source"} = [191 191 250]/255;
            col{"te_gmda"} = [205 253 254]/255;
            col{"te_wave"} = [205 253 197]/255;
            col{"sink"} = [245 194 193]/255;

            reserviors = configureDictionary("string","Box");
            [reservoirEnergy, t] = self.energyForReservoirOverTime(timeIndices=options.timeIndices);
            for iReservoir = 1:length(options.reservoirNames)
                name = options.reservoirNames(iReservoir);
                if name == "te_quadratic"
                    continue;
                end
                
                if isKey(options.customNames,name)
                    fancyName = options.customNames(name);
                else
                    fancyName = self.fancyNameForName(name);
                end

                reserviors(name) = Box(fancyName,FaceColor=col{name}, FontSize=16, CornerRadius=0.10);
                if options.shouldShowReservoirEnergy
                    energy = mean(reservoirEnergy(iReservoir).energy)/options.energy_scale;
                    flux = (reservoirEnergy(iReservoir).energy(end) - reservoirEnergy(iReservoir).energy(1))/(t(end)-t(1))/options.flux_scale;
                    if abs(flux) > options.fluxTolerance
                        if flux > 0
                            reserviors(name).Sublabel=sprintf("%.2f %s + %.2f %s",energy,options.energy_scale_units,abs(flux),options.flux_scale_units);
                        else
                            reserviors(name).Sublabel=sprintf("%.2f %s – %.2f %s",energy,options.energy_scale_units,abs(flux),options.flux_scale_units);
                        end
                    else
                        reserviors(name).Sublabel=sprintf("%.2f %s",energy,options.energy_scale_units);
                    end
                end
            end

            sources = Box.empty(0,0);
            sinks = Box.empty(0,0);
            source_arrows = Arrow.empty(0,0);
            sink_arrows = Arrow.empty(0,0);
            for iFlux=1:length(forcing_fluxes)
                if forcing_fluxes(iFlux).name == "nonlinear_advection"
                    continue;
                end
                % 
                % if abs(forcing_fluxes(iFlux).te)/options.escale < options.flux_tolerance
                %     continue;
                % end

                if isKey(options.customNames,forcing_fluxes(iFlux).name)
                    fancyName = options.customNames(forcing_fluxes(iFlux).name);
                else
                    fancyName = forcing_fluxes(iFlux).fancyName;
                end

                if forcing_fluxes(iFlux).te/options.flux_scale/2 > options.fluxTolerance
                    sources(end+1) = Box(fancyName,FaceColor=col{"source"}, FontSize=16);
                    reserviorNames = reserviors.keys;
                    for iRes=1:length(reserviorNames)
                        name = reserviorNames(iRes);
                        magnitude = abs(forcing_fluxes(iFlux).(name))/options.flux_scale;
                        if options.shouldShowUnits
                            label = sprintf("%.2f %s",magnitude,options.flux_scale_units);
                        else
                            label = sprintf("%.2f",magnitude);
                        end
                        if abs(magnitude) > options.fluxTolerance
                            source_arrows(end+1) = Arrow(sources(end),reserviors(name),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                        end
                    end
                elseif forcing_fluxes(iFlux).te/options.flux_scale/2 < -options.fluxTolerance
                    sinks(end+1) = Box(fancyName,FaceColor=col{"sink"}, FontSize=16);
                    reserviorNames = reserviors.keys;
                    for iRes=1:length(reserviorNames)
                        name = reserviorNames(iRes);
                        magnitude = abs(forcing_fluxes(iFlux).(name))/options.flux_scale;
                        if options.shouldShowUnits
                            label = sprintf("%.2f %s",magnitude,options.flux_scale_units);
                        else
                            label = sprintf("%.2f",magnitude);
                        end
                        if abs(magnitude) > options.fluxTolerance
                            sink_arrows(end+1) = Arrow(reserviors(name),sinks(end),Label=label,Magnitude=magnitude, LabelOffset=0.25, FontSize=14);
                        end
                    end
                end
            end

            % Now sort the forcing to minimize arrow crossing. The
            % reservoir order is assumed fixed. So what the heck is the
            % logic here?
            sources_sorted = Box.empty(0,0);
            reserviorNames = reserviors.keys;
            for iRes=1:length(reserviorNames)
                indices = arrayfun( @(a) a.Target == reserviors(reserviorNames(iRes)), source_arrows);
                candidate_sources = unique([source_arrows(indices).Source]);
                for k = 1:numel(candidate_sources)
                    if ~any(candidate_sources(k) == sources_sorted)
                        sources_sorted(end+1) = candidate_sources(k);  % Append if not present
                    end
                end
            end
            if length(sources) ~= length(sources_sorted)
                error("messed up my logic. this algorithm will drop sources")
            end
            sources = sources_sorted;

            inertial_arrows = Arrow.empty(0,0);
            if length(reserviors.keys) == 2
                inertial_fluxes = self.inertialFluxesForReservoirSpatialTemporalAverage(reservoirNames=options.reservoirNames,timeIndices=options.timeIndices);

                mag_geo = sum([inertial_fluxes(:).te_gmda])/options.flux_scale;
                mag_wave = sum([inertial_fluxes(:).te_wave])/options.flux_scale;
                magnitude = (abs(mag_geo) + abs(mag_wave))/2;
                if options.shouldShowUnits
                    label = sprintf("%.2f %s",magnitude,options.flux_scale_units);
                else
                    label = sprintf("%.2f",magnitude);
                end
                if magnitude > options.fluxTolerance
                    if mag_geo > 0
                        inertial_arrows(end+1) = Arrow(reserviors("te_wave"),reserviors("te_gmda"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    else
                        inertial_arrows(end+1) = Arrow(reserviors("te_gmda"),reserviors("te_wave"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    end
                end
                    
            end


            fig = plotThreeRowBoxDiagram(sources, reserviors.values, sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[3.0 1.5], Title=options.title, visible=options.visible);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (over time)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [reservoirs, t] = energyForReservoirOverTime(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te", "te_quadratic"];
                options.timeIndices = Inf;
            end
            %% Measure total energy in each reservoir
            [KE_g,PE_g,E_mda,E_w,E_io,ke,pe_quadratic,ape] =self.diagfile.readVariables('KE_g','PE_g','E_mda','E_w','E_io','ke','pe_quadratic','ape');

            if isinf(options.timeIndices)
                filter = @(v) v;
            else
                filter = @(v) v(options.timeIndices);
            end

            % we have to preallocated an array of structs
            clear reservoirs;
            reservoirs(length(options.reservoirNames)) = struct("name","placeholder");
            for iReservoir = length(options.reservoirNames):-1:1
                reservoirs(iReservoir).name = options.reservoirNames(iReservoir);
                reservoirs(iReservoir).fancyName = self.fancyNameForName(options.reservoirNames(iReservoir));
                switch options.reservoirNames(iReservoir)
                    case "ke_g"
                        reservoirs(iReservoir).energy = KE_g;
                    case "pe_g"
                        reservoirs(iReservoir).energy = PE_g;
                    case "te_g"
                        reservoirs(iReservoir).energy = KE_g + PE_g;
                    case "te_mda"
                        reservoirs(iReservoir).energy = E_mda;
                    case "te_gmda"
                        reservoirs(iReservoir).energy = KE_g + PE_g + E_mda;
                    case "te_igw"
                        reservoirs(iReservoir).energy = E_w;
                    case "te_io"
                        reservoirs(iReservoir).energy = E_io;
                    case "te_wave"
                        reservoirs(iReservoir).energy = E_w+E_io;
                    case "te"
                        reservoirs(iReservoir).energy = ke + ape;
                    case "te_quadratic"
                        reservoirs(iReservoir).energy = ke + pe_quadratic;
                    otherwise
                        error("unknown energy reservoir");
                end
                reservoirs(iReservoir).energy = filter(reservoirs(iReservoir).energy);
            end
            t = filter(self.diagfile.readVariables('t'));
        end

        function forcing_fluxes = forcingFluxesForReservoirOverTime(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
            else
                filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
            end
            forcing_fluxes = self.forcingFluxesForReservoir(reservoirNames=options.reservoirNames);
            exact_forcing_fluxes = self.exactForcingFluxesOverTime();
            for iForce=1:length(forcing_fluxes)
                forcing_fluxes(iForce).te = reshape(exact_forcing_fluxes(iForce).te,1,1,[]);
            end

            forcing_fluxes = self.filterFluxesForReservoir(forcing_fluxes,filter=filter_space);
        end

        function inertial_fluxes = inertialFluxesForReservoirOverTime(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.shouldConsolidateTriads = true;
            end
            filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxesForReservoir(reservoirNames=options.reservoirNames,shouldConsolidateTriads=options.shouldConsolidateTriads),filter=filter_space);
        end

        function forcing_fluxes = forcingFluxesForReservoirTemporalAverage(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.timeIndices = Inf;
            end

            if isinf(options.timeIndices)
                filter_space = @(v) mean(v,3);
            else
                filter_space = @(v) mean(v(:,:,options.timeIndices),3);
            end
            forcing_fluxes = self.filterFluxesForReservoir(self.forcingFluxesForReservoir(reservoirNames=options.reservoirNames),filter=filter_space);
        end

        function inertial_fluxes = inertialFluxesForReservoirTemporalAverage(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.timeIndices = Inf;
                options.shouldConsolidateTriads = true;
            end

            if isinf(options.timeIndices)
                filter_space = @(v) mean(v,3);
            else
                filter_space = @(v) mean(v(:,:,options.timeIndices),3);
            end
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxesForReservoir(reservoirNames=options.reservoirNames,shouldConsolidateTriads=options.shouldConsolidateTriads),filter=filter_space);
        end

        function forcing_fluxes = forcingFluxesForReservoirSpatialTemporalAverage(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.timeIndices = Inf;
            end

            forcing_fluxes = self.filterFluxesForReservoir(self.forcingFluxesForReservoirOverTime(reservoirNames=options.reservoirNames,timeIndices=options.timeIndices),filter=@(v) mean(v));
        end

        function inertial_fluxes = inertialFluxesForReservoirSpatialTemporalAverage(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.timeIndices = Inf;
                options.shouldConsolidateTriads = true;
            end

            if isinf(options.timeIndices)
                filter_space = @(v) sum(sum(mean(v,3),1),2);
            else
                filter_space = @(v) sum(sum(mean(v(:,:,options.timeIndices),3),1),2);
            end
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxesForReservoir(reservoirNames=options.reservoirNames,shouldConsolidateTriads=options.shouldConsolidateTriads),filter=filter_space);
        end

        function forcing_fluxes = exactForcingFluxesOverTime(self)
            arguments
                self WVDiagnostics
            end
            forcingNames = self.wvt.forcingNames;
            forcing_fluxes(length(forcingNames)) = struct("name","placeholder");

            for iForce=1:length(forcingNames)
                name = replace(forcingNames(iForce),"-","_");
                name = replace(name," ","_");
                forcing_fluxes(iForce).name = name;
                forcing_fluxes(iForce).fancyName = forcingNames(iForce);
                forcing_fluxes(iForce).te = self.diagfile.readVariables("E_" + name);
            end
        end

        function forcing_fluxes = forcingFluxesForReservoir(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
            end
            forcingNames = self.wvt.forcingNames;
            forcing_fluxes(length(forcingNames)) = struct("name","placeholder");

            for iForce=1:length(forcingNames)
                name = replace(forcingNames(iForce),"-","_");
                name = replace(name," ","_");

                % these are temporary variavbles for use within this loop only
                Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
                Ejk.Em = self.diagfile.readVariables("Em_" + name);
                Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
                Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

                forcing_fluxes(iForce).name = name;
                forcing_fluxes(iForce).fancyName = forcingNames(iForce);

                % per-reservoir fluxes
                fluxes = self.energyFluxForReservoirFromStructure(Ejk,options.reservoirNames);
                for iReservoir = 1:length(options.reservoirNames)
                    forcing_fluxes(iForce).(options.reservoirNames(iReservoir)) = fluxes{iReservoir};
                end
            end
        end

        function inertial_fluxes = inertialFluxesForReservoir(self,options)
            arguments
                self WVDiagnostics
                options.reservoirNames = ["te_gmda", "te_wave", "te_quadratic"];
                options.shouldConsolidateTriads = true;
            end

            triadFlowComponents_t = self.wvt.primaryFlowComponents;
            inertial_fluxes(length(triadFlowComponents_t)^2) = struct("name","placeholder");

            counter = 1;
            for i=1:length(triadFlowComponents_t)
                for j=1:length(triadFlowComponents_t)
                    name = triadFlowComponents_t(i).abbreviatedName + "_" + triadFlowComponents_t(j).abbreviatedName;

                    % these are temporary variavbles for use within this loop only
                    Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
                    Ejk.Em = self.diagfile.readVariables("Em_" + name);
                    Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
                    Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

                    % total fluxes
                    inertial_fluxes(counter).name = name;
                    inertial_fluxes(counter).fancyName = triadFlowComponents_t(i).shortName + "-" + triadFlowComponents_t(j).shortName;

                    % per-reservoir fluxes
                    fluxes = self.energyFluxForReservoirFromStructure(Ejk,options.reservoirNames);
                    for iReservoir = 1:length(options.reservoirNames)
                        inertial_fluxes(counter).(options.reservoirNames(iReservoir)) = fluxes{iReservoir};
                    end

                    % increment counter
                    counter = counter+1;
                end
            end

            if options.shouldConsolidateTriads
                arraySize = size(inertial_fluxes(counter-1).(options.reservoirNames(iReservoir)));

                consolidatedFlowComponents_t = WVFlowComponent.empty(0,0);
                for iReservoir = 1:length(options.reservoirNames)
                    switch options.reservoirNames(iReservoir)
                        case "ke_g"
                            consolidatedFlowComponents_t(end+1) = self.wvt.geostrophicComponent;
                        case "pe_g"
                            consolidatedFlowComponents_t(end+1) = self.wvt.geostrophicComponent;
                        case "te_g"
                            consolidatedFlowComponents_t(end+1) = self.wvt.geostrophicComponent;
                        case "te_mda"
                            consolidatedFlowComponents_t(end+1) = self.wvt.mdaComponent;
                        case "te_gmda"
                            consolidatedFlowComponents_t(end+1) = self.wvt.geostrophicComponent + self.wvt.mdaComponent;
                        case "te_igw"
                            consolidatedFlowComponents_t(end+1) = self.wvt.waveComponent;
                        case "te_io"
                            consolidatedFlowComponents_t(end+1) = self.wvt.inertialComponent;
                        case "te_wave"
                            consolidatedFlowComponents_t(end+1) = self.wvt.waveComponent + self.wvt.inertialComponent;
                    end
                end
                % consolidatedFlowComponents_t = unique(consolidatedFlowComponents_t);
                inertial_fluxes_consol(length(consolidatedFlowComponents_t)^2) = struct("name","placeholder");
                n = length(consolidatedFlowComponents_t);
                for i=1:length(consolidatedFlowComponents_t)
                    for j=1:length(consolidatedFlowComponents_t)
                        counter = j+(i-1)*n;
                        inertial_fluxes_consol(counter).name = consolidatedFlowComponents_t(i).abbreviatedName + "_" + consolidatedFlowComponents_t(j).abbreviatedName;
                        inertial_fluxes_consol(counter).fancyName = consolidatedFlowComponents_t(i).shortName + "-" + consolidatedFlowComponents_t(j).shortName;
                        for iReservoir = 1:length(options.reservoirNames)
                            inertial_fluxes_consol(counter).(options.reservoirNames(iReservoir)) = zeros(arraySize);
                        end
                    end
                end

                counter = 1;
                for i=1:length(triadFlowComponents_t)
                    index_i = find(arrayfun(@(a) a.contains(triadFlowComponents_t(i)), consolidatedFlowComponents_t));
                    for j=1:length(triadFlowComponents_t)
                        index_j = find(arrayfun(@(a) a.contains(triadFlowComponents_t(j)), consolidatedFlowComponents_t));
                        for iReservoir = 1:length(options.reservoirNames)
                            inertial_fluxes_consol(index_j+(index_i-1)*n).(options.reservoirNames(iReservoir)) = inertial_fluxes(counter).(options.reservoirNames(iReservoir)) + inertial_fluxes_consol(index_j+(index_i-1)*n).(options.reservoirNames(iReservoir));
                        end
                        counter = counter+1;
                    end
                end

                inertial_fluxes = inertial_fluxes_consol;
            end
        end
    end

    methods (Static)
        cmap = cmocean(ColormapName,varargin)

        function iTimeChanged(metaProp,eventData)
            wvd = eventData.AffectedObject;
            wvd.wvt.initFromNetCDFFile(wvd.wvfile,iTime=wvd.iTime);
        end

        function fancyName = fancyNameForName(name)
            switch name
                case "ke_g"
                    fancyName = "geostrophic kinetic";
                case "pe_g"
                    fancyName = "geostrophic potential";
                case "te_g"
                    fancyName = "geostrophic";
                case "te_mda"
                    fancyName = "mean density anomaly";
                case "te_gmda"
                    fancyName = "geostrophic + mda";
                case "te_igw"
                    fancyName = "internal gravity wave";
                case "te_io"
                    fancyName = "inertial";
                case "te_wave"
                    fancyName = "wave";
                case "te"
                    fancyName = "total";
                case "te_quadratic"
                    fancyName = "total quadratic";
                otherwise
                    error("unknown energy reservoir");
            end
        end

        function filtered_fluxes = filterFluxesForReservoir(fluxes,options)
            arguments
                fluxes
                options.filter = @(v) v;
            end

            reservoirNames = string(fields(fluxes));
            reservoirNames(reservoirNames=="name") = [];
            reservoirNames(reservoirNames=="fancyName") = [];
            for iForce=1:length(fluxes)
                for iReservoir = 1:length(reservoirNames)
                    fluxes(iForce).(reservoirNames(iReservoir)) = options.filter(fluxes(iForce).(reservoirNames(iReservoir)));
                end
            end
            filtered_fluxes = fluxes;
        end

        function eFlux = energyFluxForReservoirFromStructure(Ejk,reservoirNames)
            eFlux = cell(length(reservoirNames),1);
            for iReservoir = 1:length(reservoirNames)
                switch reservoirNames(iReservoir)
                    case "ke_g"
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:);
                    case "pe_g"
                        eFlux{iReservoir} = Ejk.PE0(:,2:end,:);
                    case "te_g"
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:)+Ejk.PE0(:,2:end,:);
                    case "te_mda"
                        eFlux{iReservoir} = Ejk.PE0(:,1,:);
                    case "te_gmda"
                        eFlux{iReservoir} = Ejk.KE0 + Ejk.PE0;
                    case "te_igw"
                        eFlux{iReservoir} = Ejk.Ep(:,2:end,:) + Ejk.Em(:,2:end,:);
                    case "te_io"
                        eFlux{iReservoir} = Ejk.Ep(:,1,:) + Ejk.Em(:,1,:);
                    case "te_wave"
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em;
                    case "te_quadratic"
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em+Ejk.KE0+Ejk.PE0;
                    otherwise
                        error("unknown energy reservoir");
                end
            end
        end
    end
end