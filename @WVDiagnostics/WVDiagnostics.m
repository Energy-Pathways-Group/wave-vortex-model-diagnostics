classdef WVDiagnostics < handle
    %WVDiagnostics Produces diagnostics and figures from WVModel output
    %   This is a collection of diagnostic tools for analyzing model output
    properties
        wvpath
        diagpath
        wvfile
        diagfile
        wvt
    end

    properties (Dependent)
        t_diag
        t_wv
        j
        kRadial
    end

    properties (SetObservable)
        iTime = Inf
    end

    methods
        function self = WVDiagnostics(filename,options)
            % Initializes the WVDiagnostics object, loads the wave-vortex transform and diagnostics files.
            %
            % - Topic: Constructor
            % - Declaration: self = WVDiagnostics(filename,options)
            % - Parameter filename: path to the WVModel output file
            % - Parameter options.diagnosticsFilePath: (optional) path to the diagnostics file
            % - Returns self: the constructed WVDiagnostics object
            arguments
                filename
                options.diagnosticsFilePath
            end
            self.wvpath = filename;


            [self.wvt, self.wvfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=Inf);
            [fpath,fname,~] = fileparts(filename);
            if ~isfield(options,"diagnosticsFilePath")
                if ~isempty(fpath)
                    self.diagpath = fullfile(fpath,strcat(fname,"-diagnostics.nc"));
                else
                    self.diagpath= fullfile(strcat(fname,"-diagnostics.nc"));
                end
            else
                self.diagpath = options.diagnosticsFilePath;
            end
            if exist(self.diagpath,"file")
                self.diagfile = NetCDFFile(self.diagpath);
            else
                warning("No diagnostics file found. Some functionality will not be available.")
            end

            addlistener(self,'iTime','PostSet',@WVDiagnostics.iTimeChanged);
        end

        function t = get.t_diag(self)
            % Get time vector from the diagnostics file
            %
            % Reads the 't' variable from the diagnostics file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.t_diag(self)
            % - Returns t: time vector from diagnostics file
            t = self.diagfile.readVariables('t');
        end

        function t = get.j(self)
            % Get j indices
            %
            % Reads the 'j' variable from the diagnostics file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.j(self)
            % - Returns t: j indices from diagnostics file
            t = self.diagfile.readVariables('j');
        end

        function t = get.kRadial(self)
            % Get kRadial indices
            %
            % Reads the 'kRadial' variable from the diagnostics file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.kRadial(self)
            % - Returns t: kRadial indices from diagnostics file
            t = self.diagfile.readVariables('kRadial');
        end

        function t = get.t_wv(self)
            % Get time vector from the model output file
            %
            % Reads the 't' variable from the wave-vortex file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.t_wv(self)
            % - Returns t: time vector from wave-vortex file
            t = self.wvfile.readVariables('t');
        end

        createDiagnosticsFile(self,options)
        fig = plotMooringRotarySpectrum(self)
        fig = plotFluidStateMultipanel(self,options)

        function fig = plotEnstrophyOverTime(self,options)
            % Plot enstrophy over time
            %
            % Plots the quadratic and APV enstrophy as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotEnstrophyOverTime(self,options)
            % - Parameter options.tscale: time scaling factor (default: 86400)
            % - Parameter options.tscale_units: units for time axis (default: "days")
            % - Parameter options.visible: figure visibility (default: "on")
            % - Returns fig: handle to the generated figure
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
            % Plot energy for each reservoir over time
            %
            % Plots the energy in each specified reservoir as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotEnergyForReservoirOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.shouldIncludeExactTotalEnergy: include exact total energy (default: true)
            % - Parameter options.tscale: time scaling factor (default: 86400)
            % - Parameter options.tscale_units: units for time axis (default: "days")
            % - Parameter options.escale: energy scaling factor (default: 3.74)
            % - Parameter options.escale_units: units for energy axis (default: "GM")
            % - Parameter options.visible: figure visibility (default: "on")
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.shouldIncludeExactTotalEnergy = true
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74
                options.escale_units = "GM"
                options.visible = "on"
            end
            reservoirs = self.energyForReservoirOverTime(energyReservoirs=options.energyReservoirs,shouldIncludeExactTotalEnergy=options.shouldIncludeExactTotalEnergy);
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
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxForReservoirOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.tscale: time scaling factor (default: 86400)
            % - Parameter options.tscale_units: units for time axis (default: "days")
            % - Parameter options.escale: energy scaling factor (default: 3.74/(86400*365))
            % - Parameter options.escale_units: units for energy axis (default: "GM/yr")
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74/(86400*365)
                options.escale_units = "GM/yr"
                options.visible = "on"
                options.filter = @(v) v;
            end
            forcing_fluxes = self.forcingFluxesOverTime(energyReservoirs=options.energyReservoirs);
            t = self.diagfile.readVariables('t');

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
            for iReservoir = 1:length(options.energyReservoirs)
                nexttile(tl);
                for iForce = 1:length(forcing_fluxes)
                    plot(t/options.tscale,options.filter(forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/options.escale)), hold on
                end
                legend(forcing_fluxes.fancyName)

                fancyName = options.energyReservoirs(iReservoir).fancyName;
                xlabel("time (" + options.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + options.escale_units + ")")
                xlim([min(t) max(t)]/options.tscale);
            end
        end

        function fig = plotInertialFluxForReservoirOverTime(self,options)
            % Plot inertial flux for each reservoir over time
            %
            % Plots the energy flux between reservoirs due to inertial interactions as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotInertialFluxForReservoirOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.tscale: time scaling factor (default: 86400)
            % - Parameter options.tscale_units: units for time axis (default: "days")
            % - Parameter options.escale: energy scaling factor (default: 3.74/(86400*365))
            % - Parameter options.escale_units: units for energy axis (default: "GM/yr")
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.tscale = 86400
                options.tscale_units = "days"
                options.escale = 3.74/(86400*365)
                options.escale_units = "GM/yr"
                options.visible = "on"
                options.filter = @(v) v;
            end
            inertial_fluxes = self.inertialFluxesOverTime(energyReservoirs=options.energyReservoirs);
            t = self.diagfile.readVariables('t');

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
            for iReservoir = 1:length(options.energyReservoirs)
                nexttile(tl);
                for iForce = 1:length(inertial_fluxes)
                    plot(t/options.tscale,options.filter(inertial_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/options.escale)), hold on
                end
                legend(inertial_fluxes.fancyName)

                fancyName = options.energyReservoirs(iReservoir).fancyName;
                xlabel("time (" + options.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + options.escale_units + ")")
                xlim([min(t) max(t)]/options.tscale);
            end
        end

        function fig = plotSourcesSinksReservoirsDiagram(self,options)
            % Plot sources, sinks, and reservoirs diagram
            %
            % Generates a diagram showing energy sources, sinks, and reservoirs, including fluxes between them.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotSourcesSinksReservoirsDiagram(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave])
            % - Parameter options.flux_scale: scaling for fluxes (default: 3.74/(86400*365))
            % - Parameter options.flux_scale_units: units for fluxes (default: "GM/yr")
            % - Parameter options.energy_scale: scaling for energy (default: 3.74)
            % - Parameter options.energy_scale_units: units for energy (default: "GM")
            % - Parameter options.customNames: dictionary for custom names
            % - Parameter options.fluxTolerance: tolerance for displaying fluxes (default: 1e-2)
            % - Parameter options.shouldShowUnits: show units in labels (default: true)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Parameter options.shouldShowReservoirEnergy: show reservoir energy (default: true)
            % - Parameter options.title: diagram title (default: "Energy Pathways")
            % - Parameter options.visible: figure visibility (default: "on")
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave];
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
            forcing_fluxes = self.forcingFluxesSpatialTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

            col = configureDictionary("string","cell");
            col{"source"} = [191 191 250]/255;
            col{"ke_g"} = [205 253 254]/255;
            col{"pe_g"} = [205 253 254]/255;
            col{"te_gmda"} = [205 253 254]/255;
            col{"te_wave"} = [205 253 197]/255;
            col{"sink"} = [245 194 193]/255;

            reservoirs = configureDictionary("string","Box");
            [reservoirEnergy, t] = self.energyForReservoirOverTime(timeIndices=options.timeIndices);
            for iReservoir = 1:length(options.energyReservoirs)
                name = options.energyReservoirs(iReservoir).name;
                if name == "te_quadratic"
                    continue;
                end
                
                if isKey(options.customNames,name)
                    fancyName = options.customNames(name);
                else
                    fancyName = self.fancyNameForName(name);
                end

                reservoirs(name) = Box(fancyName,FaceColor=col{name}, FontSize=16, CornerRadius=0.10);
                if options.shouldShowReservoirEnergy
                    energy = mean(reservoirEnergy(iReservoir).energy)/options.energy_scale;
                    flux = (reservoirEnergy(iReservoir).energy(end) - reservoirEnergy(iReservoir).energy(1))/(t(end)-t(1))/options.flux_scale;
                    if abs(flux) > options.fluxTolerance
                        if flux > 0
                            reservoirs(name).Sublabel=sprintf("%.2f %s + %.2f %s",energy,options.energy_scale_units,abs(flux),options.flux_scale_units);
                        else
                            reservoirs(name).Sublabel=sprintf("%.2f %s – %.2f %s",energy,options.energy_scale_units,abs(flux),options.flux_scale_units);
                        end
                    else
                        reservoirs(name).Sublabel=sprintf("%.2f %s",energy,options.energy_scale_units);
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
                    reservoirNames = reservoirs.keys;
                    for iRes=1:length(reservoirNames)
                        name = reservoirNames(iRes);
                        magnitude = abs(forcing_fluxes(iFlux).(name))/options.flux_scale;
                        if options.shouldShowUnits
                            label = sprintf("%.2f %s",magnitude,options.flux_scale_units);
                        else
                            label = sprintf("%.2f",magnitude);
                        end
                        if abs(magnitude) > options.fluxTolerance
                            source_arrows(end+1) = Arrow(sources(end),reservoirs(name),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                        end
                    end
                elseif forcing_fluxes(iFlux).te/options.flux_scale/2 < -options.fluxTolerance
                    sinks(end+1) = Box(fancyName,FaceColor=col{"sink"}, FontSize=16);
                    reservoirNames = reservoirs.keys;
                    for iRes=1:length(reservoirNames)
                        name = reservoirNames(iRes);
                        magnitude = abs(forcing_fluxes(iFlux).(name))/options.flux_scale;
                        if options.shouldShowUnits
                            label = sprintf("%.2f %s",magnitude,options.flux_scale_units);
                        else
                            label = sprintf("%.2f",magnitude);
                        end
                        if abs(magnitude) > options.fluxTolerance
                            sink_arrows(end+1) = Arrow(reservoirs(name),sinks(end),Label=label,Magnitude=magnitude, LabelOffset=0.25, FontSize=14);
                        end
                    end
                end
            end

            % Now sort the forcing to minimize arrow crossing. The
            % reservoir order is assumed fixed. So what the heck is the
            % logic here?
            sources_sorted = Box.empty(0,0);
            reservoirNames = reservoirs.keys;
            for iRes=1:length(reservoirNames)
                indices = arrayfun( @(a) a.Target == reservoirs(reservoirNames(iRes)), source_arrows);
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
            if length(reservoirs.keys) == 2
                inertial_fluxes = self.inertialFluxesSpatialTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

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
                        inertial_arrows(end+1) = Arrow(reservoirs("te_wave"),reservoirs("te_gmda"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    else
                        inertial_arrows(end+1) = Arrow(reservoirs("te_gmda"),reservoirs("te_wave"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    end
                end
                    
            end


            fig = plotThreeRowBoxDiagram(sources, reservoirs.values, sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[3.0 1.5], Title=options.title, visible=options.visible);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (over time)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [reservoirs, t] = energyForReservoirOverTime(self,options)
            % Compute energy for each reservoir over time
            %
            % Returns the energy in each specified reservoir as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: [reservoirs, t] = energyForReservoirOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.shouldIncludeExactTotalEnergy: include exact total energy (default: true)
            % - Parameter options.timeIndices: indices for time selection (default: Inf)
            % - Returns reservoirs: struct array with energy for each reservoir
            % - Returns t: time vector
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.shouldIncludeExactTotalEnergy = true
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
            if options.shouldIncludeExactTotalEnergy
                reservoirs(length(options.energyReservoirs)+1) = struct("name","placeholder");
            else
                reservoirs(length(options.energyReservoirs)) = struct("name","placeholder");
            end
            for iReservoir = 1:length(options.energyReservoirs)
                reservoirs(iReservoir).name = options.energyReservoirs(iReservoir).name;
                reservoirs(iReservoir).fancyName = options.energyReservoirs(iReservoir).fancyName;
                switch options.energyReservoirs(iReservoir)
                    case EnergyReservoir.geostrophic_kinetic
                        reservoirs(iReservoir).energy = KE_g;
                    case EnergyReservoir.geostrophic_potential
                        reservoirs(iReservoir).energy = PE_g;
                    case EnergyReservoir.geostrophic
                        reservoirs(iReservoir).energy = KE_g + PE_g;
                    case EnergyReservoir.mda
                        reservoirs(iReservoir).energy = E_mda;
                    case EnergyReservoir.geostrophic_mda
                        reservoirs(iReservoir).energy = KE_g + PE_g + E_mda;
                    case EnergyReservoir.igw
                        reservoirs(iReservoir).energy = E_w;
                    case EnergyReservoir.io
                        reservoirs(iReservoir).energy = E_io;
                    case EnergyReservoir.wave
                        reservoirs(iReservoir).energy = E_w+E_io;
                    case EnergyReservoir.total
                        reservoirs(iReservoir).energy = ke + pe_quadratic;
                    otherwise
                        error("unknown energy reservoir");
                end
                reservoirs(iReservoir).energy = filter(reservoirs(iReservoir).energy);
            end
            if options.shouldIncludeExactTotalEnergy
                reservoirs(end).name = "te";
                reservoirs(end).fancyName = "total";
                reservoirs(end).energy = ke + ape;
                reservoirs(end).energy = filter(reservoirs(end).energy);
            end

            t = filter(self.diagfile.readVariables('t'));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Flux averages, scalar
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function forcing_fluxes = forcingFluxesSpatialTemporalAverage(self,options)
            % Compute spatial-temporal average of forcing fluxes
            %
            % Returns the spatial-temporal average of energy fluxes from external forcing for each reservoir.
            %
            % - Topic: Flux averages, scalar
            % - Declaration: forcing_fluxes = forcingFluxesSpatialTemporalAverage(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns forcing_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
            end

            forcing_fluxes = self.filterFluxesForReservoir(self.forcingFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices),filter=@(v) mean(v));
        end

        function inertial_fluxes = inertialFluxesSpatialTemporalAverage(self,options)
            % Compute spatial-temporal average of inertial fluxes
            %
            % Returns the spatial-temporal average of energy fluxes due to inertial interactions for each reservoir.
            %
            % - Topic: Flux averages, scalar
            % - Declaration: inertial_fluxes = inertialFluxesSpatialTemporalAverage(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Parameter options.triadComponents: vector of TriadFlowComponent objects (default: [geostrophic_mda, wave])
            % - Returns inertial_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
            end

            if isinf(options.timeIndices)
                filter_space = @(v) sum(sum(mean(v,3),1),2);
            else
                filter_space = @(v) sum(sum(mean(v(:,:,options.timeIndices),3),1),2);
            end
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes over time, [t 1]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function forcing_fluxes = exactForcingFluxesOverTime(self)
            % Compute exact forcing fluxes over time
            %
            % Returns the exact energy fluxes from external forcing for each time step.
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: forcing_fluxes = exactForcingFluxesOverTime(self)
            % - Returns forcing_fluxes: struct array with exact fluxes
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

        function forcing_fluxes = forcingFluxesOverTime(self,options)
            % Compute forcing fluxes over time
            %
            % Returns the energy fluxes from external forcing for each reservoir as a function of time.
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: forcing_fluxes = forcingFluxesOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.timeIndices: indices for time selection (default: Inf)
            % - Returns forcing_fluxes: struct array with fluxes over time
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
            else
                filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
            end
            forcing_fluxes = self.forcingFluxes(energyReservoirs=options.energyReservoirs);
            exact_forcing_fluxes = self.exactForcingFluxesOverTime();
            for iForce=1:length(forcing_fluxes)
                forcing_fluxes(iForce).te = reshape(exact_forcing_fluxes(iForce).te,1,1,[]);
            end

            forcing_fluxes = self.filterFluxesForReservoir(forcing_fluxes,filter=filter_space);
        end

        function inertial_fluxes = inertialFluxesOverTime(self,options)
            % Compute inertial fluxes over time
            %
            % Returns the energy fluxes due to inertial interactions for each reservoir as a function of time.
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: inertial_fluxes = inertialFluxesOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.triadComponents: vector of TriadFlowComponent objects (default: [geostrophic_mda, wave])
            % - Returns inertial_fluxes: struct array with fluxes over time
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
            end
            filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes in space, [j kRadial]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function forcing_fluxes = forcingFluxesTemporalAverage(self,options)
            % Compute temporally averaged forcing fluxes
            %
            % Returns the temporally averaged energy fluxes from external forcing for each reservoir.
            %
            % - Topic: Fluxes in space, [j kRadial]
            % - Declaration: forcing_fluxes = forcingFluxesTemporalAverage(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns forcing_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
            end

            if isinf(options.timeIndices)
                filter_space = @(v) mean(v,3);
            else
                filter_space = @(v) mean(v(:,:,options.timeIndices),3);
            end
            forcing_fluxes = self.filterFluxesForReservoir(self.forcingFluxes(energyReservoirs=options.energyReservoirs),filter=filter_space);
        end


        function inertial_fluxes = inertialFluxesTemporalAverage(self,options)
            % Computes the temporally averaged inertial fluxes.
            %
            % Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial].
            %
            % - Topic: Fluxes in space, [j kRadial]
            % - Declaration: inertial_fluxes = inertialFluxesTemporalAverage(options)
            % - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
            % - Parameter timeIndices: (optional) indices specifying which time steps to average over. Defaults to Inf (all).
            % - Parameter triadComponents: (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave].
            % - Returns inertial_fluxes: an array of structs
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
            end

            if isinf(options.timeIndices)
                filter_space = @(v) mean(v,3);
            else
                filter_space = @(v) mean(v(:,:,options.timeIndices),3);
            end
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Core fluxes functions, [j kRadial t]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function forcing_fluxes = forcingFluxes(self,options)
            % Return the energy flux from the forcing terms
            %
            % Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
            %
            % - Topic: Core function — spatial temporal
            % - Declaration: forcing_fluxes = forcingFluxes(options)
            % - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
            % - Returns forcing_fluxes: an array of structs
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
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
                fluxes = self.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
                for iReservoir = 1:length(options.energyReservoirs)
                    forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name) = fluxes{iReservoir};
                end
            end
        end

        function inertial_fluxes = inertialFluxes(self,options)
            % Return the energy flux from the inertial terms, specified as triad components
            %
            % Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
            %
            % - Topic: Core function — spatial temporal
            % - Declaration: inertial_fluxes = inertialFluxes(options)
            % - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
            % - Parameter triadComponents: (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave].
            % - Returns inertial_fluxes: an array of structs
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
            end

            primaryFlowComponents_t = self.wvt.primaryFlowComponents;
            inertial_fluxes(length(primaryFlowComponents_t)^2) = struct("name","placeholder");

            counter = 1;
            for i=1:length(primaryFlowComponents_t)
                for j=1:length(primaryFlowComponents_t)
                    name = primaryFlowComponents_t(i).abbreviatedName + "_" + primaryFlowComponents_t(j).abbreviatedName;

                    % these are temporary variavbles for use within this loop only
                    Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
                    Ejk.Em = self.diagfile.readVariables("Em_" + name);
                    Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
                    Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

                    % total fluxes
                    inertial_fluxes(counter).name = name;
                    inertial_fluxes(counter).fancyName = primaryFlowComponents_t(i).shortName + "-" + primaryFlowComponents_t(j).shortName;

                    % per-reservoir fluxes
                    fluxes = self.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
                    for iReservoir = 1:length(options.energyReservoirs)
                        inertial_fluxes(counter).(options.energyReservoirs(iReservoir).name) = fluxes{iReservoir};
                    end

                    % increment counter
                    counter = counter+1;
                end
            end

            consolidatedFlowComponents_t = WVFlowComponent.empty(0,0);
            for i=1:length(options.triadComponents)
                consolidatedFlowComponents_t(end+1) = options.triadComponents(i).flowComponent(self.wvt);
            end

            inertial_fluxes_consol(length(options.triadComponents)^2) = struct("name","placeholder");
            n = length(options.triadComponents);
            for i=1:n
                for j=1:n
                    counter = j+(i-1)*n;
                    inertial_fluxes_consol(counter).name = options.triadComponents(i).name + "_" + options.triadComponents(j).name;
                    inertial_fluxes_consol(counter).fancyName = options.triadComponents(i).name + "{\nabla}" + options.triadComponents(j).name;
                    for iReservoir = 1:length(options.energyReservoirs)
                        arraySize = size(inertial_fluxes(1).(options.energyReservoirs(iReservoir).name));
                        inertial_fluxes_consol(counter).(options.energyReservoirs(iReservoir).name) = zeros(arraySize);
                    end
                end
            end

            counter = 1;
            for i=1:length(primaryFlowComponents_t)
                index_i = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(i)), consolidatedFlowComponents_t));
                for j=1:length(primaryFlowComponents_t)
                    index_j = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(j)), consolidatedFlowComponents_t));
                    for iReservoir = 1:length(options.energyReservoirs)
                        a = inertial_fluxes(counter).(options.energyReservoirs(iReservoir).name);
                        b = inertial_fluxes_consol(index_j+(index_i-1)*n).(options.energyReservoirs(iReservoir).name);
                        inertial_fluxes_consol(index_j+(index_i-1)*n).(options.energyReservoirs(iReservoir).name) = a + b;
                    end
                    counter = counter+1;
                end
            end

            inertial_fluxes = inertial_fluxes_consol;

        end
    end

    methods (Static)
        cmap = cmocean(ColormapName,varargin)

        function iTimeChanged(metaProp,eventData)
            wvd = eventData.AffectedObject;
            wvd.wvt.initFromNetCDFFile(wvd.wvfile,iTime=wvd.iTime);
        end

        function bool = areEnergyReservoirsComplete(reservoirs)
            bool = all(sum([reservoirs.vectorContents],2)==1);
            if ~bool
                warning('The collection of energy reservoirs is either not complete or over complete')
            end
        end

        function bool = areTriadComponentsComplete(reservoirs)
            bool = all(sum([reservoirs.vectorContents],2)==1);
            if ~bool
                warning('The collection of triad components is either not complete or over complete')
            end
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
                    case EnergyReservoir.geostrophic_kinetic
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:);
                    case EnergyReservoir.geostrophic_potential
                        eFlux{iReservoir} = Ejk.PE0(:,2:end,:);
                    case EnergyReservoir.geostrophic
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:)+Ejk.PE0(:,2:end,:);
                    case EnergyReservoir.mda
                        eFlux{iReservoir} = Ejk.PE0(:,1,:);
                    case EnergyReservoir.geostrophic_mda
                        eFlux{iReservoir} = Ejk.KE0 + Ejk.PE0;
                    case EnergyReservoir.igw
                        eFlux{iReservoir} = Ejk.Ep(:,2:end,:) + Ejk.Em(:,2:end,:);
                    case EnergyReservoir.io
                        eFlux{iReservoir} = Ejk.Ep(:,1,:) + Ejk.Em(:,1,:);
                    case EnergyReservoir.wave
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em;
                    case EnergyReservoir.total
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em+Ejk.KE0+Ejk.PE0;
                    otherwise
                        error("unknown energy reservoir");
                end
            end
        end
    end
end