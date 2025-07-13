classdef WVDiagnostics < handle
    %WVDiagnostics Produces diagnostics and figures from WVModel output
    %   This is a collection of diagnostic tools for analyzing model output
    properties
        wvpath
        diagpath
        wvfile
        diagfile
        wvt

        % Default scaling and units for time, energy, and flux
        tscale = 86400
        tscale_units = "days"
        escale = 3.74
        escale_units = "GM"
        flux_scale = 3.74/(86400*365)
        flux_scale_units = "GM/yr"
        z_flux_scale
        z_flux_scale_units = "m f^2/yr"
        zscale
        zscale_units = "m f^2";
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
            self.wvt.addOperation(EtaTrueOperation());
            self.wvt.addOperation(APEOperation(self.wvt));
            self.wvt.addOperation(APVOperation());
            self.wvt.addOperation(SpatialForcingOperation(self.wvt));

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
            self.zscale = self.wvt.f^2;
            self.z_flux_scale = self.zscale/(86400*365);

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
            t = self.wvfile.readVariables('wave-vortex/t');
        end

        function setEnergyUnits(self, units)
            % Set the time and energy scaling and units for plotting and output.
            %
            % Sets tscale, tscale_units, escale, and escale_units based on the specified units.
            %
            % - Topic: Configuration
            % - Declaration: setEnergyUnits(self, units)
            % - Parameter units: must be one of "si", "gm", or "si-yr"
            % - Returns: None
            arguments
                self
                units {mustBeMember(units, ["si", "gm", "si-yr"])}
            end
            switch lower(units)
                case "si"
                    self.escale = 1;
                    self.escale_units = "m^3 s^{-2}";
                    self.flux_scale = 1;
                    self.flux_scale_units = "m^3 s^{-3}";
                case "gm"
                    self.escale = 3.74;
                    self.escale_units = "GM";
                    self.flux_scale = 3.74/(86400*365);
                    self.flux_scale_units = "GM/yr";
                case "si-yr"
                    self.escale = 1;
                    self.escale_units = "m^3 s^{-2}";
                    self.flux_scale = 1/(86400*365);
                    self.flux_scale_units = "m^3 s^{-2} yr^{-1}";
                otherwise
                    error("Unknown units: %s. Must be 'si', 'gm', or 'si-yr'.", units);
            end
        end

        createDiagnosticsFile(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (instantaneous snapshot, from model output)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotFluidStateMultipanel(self,options)
        fig = plotFluidDecompositionMultipanel(self,options)
        fig = plotEnstrophySpectrum(self,options)
        fig = plotEnergySpectrum(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (ancillary data, from model output)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotMooringRotarySpectrum(self)

        fig = plotEnergyFluxTemporalAverage(self,options)
        
        function fig = plotEnstrophyOverTime(self,options)
            % Plot enstrophy over time
            %
            % Plots the quadratic and APV enstrophy as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotEnstrophyOverTime(self,options)
            % - Parameter options.visible: figure visibility (default: "on")
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.visible = "on"
                options.timeIndices = Inf;
            end
            [Z_quadratic, t] = self.enstrophyQGPVOverTime(timeIndices=options.timeIndices);
            [Z_apv, ~] = self.enstrophyAPVOverTime(timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            plot(t/self.tscale,Z_quadratic/self.zscale,LineWidth=2), hold on
            plot(t/self.tscale,Z_apv/self.zscale,LineWidth=2)
            legend('quadratic','apv')

            xlabel("time (" + self.tscale_units + ")")
            ylabel("enstrophy (" + self.zscale_units + ")")
            xlim([min(t) max(t)]/self.tscale);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (over time, from diagnostics file)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function fig = plotEnergyOverTime(self,options)
            % Plot energy for each reservoir over time
            %
            % Plots the energy in each specified reservoir as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotEnergyOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.shouldIncludeExactTotalEnergy: include exact total energy (default: true)
            % - Parameter options.visible: figure visibility (default: "on")
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.shouldIncludeExactTotalEnergy = true
                options.timeIndices = Inf;
                options.visible = "on"
            end
            [reservoirs, t] = self.energyOverTime(energyReservoirs=options.energyReservoirs,shouldIncludeExactTotalEnergy=options.shouldIncludeExactTotalEnergy,timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);

            for iReservoir = 1:length(reservoirs)
                switch reservoirs(iReservoir).name
                    case "te"
                        plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2, Color=[0 0 0]), hold on
                    case "te_quadratic"
                        plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2, Color=[0 0 0], LineStyle="-."), hold on
                    otherwise
                        plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2), hold on
                end
            end
            legend(reservoirs.fancyName)

            xlabel("time (" + self.tscale_units + ")")
            ylabel("energy (" + self.escale_units + ")")
            xlim([min(t) max(t)]/self.tscale);
        end

        function fig = plotForcingFluxOverTime(self,options)
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v) v;
            end
            [forcing_fluxes, t] = self.forcingFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
            for iReservoir = 1:length(options.energyReservoirs)
                nexttile(tl);
                for iForce = 1:length(forcing_fluxes)
                    plot(t/self.tscale,options.filter(forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/self.flux_scale)), hold on
                end
                legend(forcing_fluxes.fancyName)

                fancyName = options.energyReservoirs(iReservoir).fancyName;
                xlabel("time (" + self.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + self.flux_scale_units + ")")
                xlim([min(t) max(t)]/self.tscale);
            end
        end

        function fig = plotExactForcingFluxOverTime(self,options)
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v) v;
            end
            [forcing_fluxes, t] = self.exactForcingFluxesOverTime(timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(1,1,TileSpacing="compact");
            for iForce = 1:length(forcing_fluxes)
                plot(t/self.tscale,options.filter(forcing_fluxes(iForce).te/self.flux_scale)), hold on
            end
            legend(forcing_fluxes.fancyName)

            xlabel("time (" + self.tscale_units + ")")
            ylabel("flux (" + self.flux_scale_units + ")")
            xlim([min(t) max(t)]/self.tscale);
        end

        function fig = plotInertialFluxOverTime(self,options)
            % Plot inertial flux for each reservoir over time
            %
            % Plots the energy flux between reservoirs due to inertial interactions as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotInertialFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave];
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v) v;
            end
            [inertial_fluxes,t] = self.inertialFluxesOverTime(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
            for iReservoir = 1:length(options.energyReservoirs)
                nexttile(tl);
                for iForce = 1:length(inertial_fluxes)
                    plot(t/self.tscale,options.filter(inertial_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/self.flux_scale)), hold on
                end
                legend(inertial_fluxes.fancyName)

                fancyName = options.energyReservoirs(iReservoir).fancyName;
                xlabel("time (" + self.tscale_units + ")")
                ylabel("flux into " + fancyName + " (" + self.flux_scale_units + ")")
                xlim([min(t) max(t)]/self.tscale);
            end
        end

        function fig = plotExactEnstrophyFluxOverTime(self,options)
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v) v;
            end
            [forcing_fluxes, t] = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(1,1,TileSpacing="compact");
            for iForce = 1:length(forcing_fluxes)
                plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale)), hold on
            end
            legend(forcing_fluxes.fancyName)

            xlabel("time (" + self.tscale_units + ")")
            ylabel("flux (" + self.z_flux_scale_units + ")")
            xlim([min(t) max(t)]/self.tscale);
        end

        function fig = plotEnstrophyFluxOverTime(self,options)
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % Some good filters are:
            % filter=@(v,t) movmean(v,21);
            % filter=@(v,t) cumtrapz(t,v)./(t+1)
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v,t) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v,t) v;
            end
            [forcing_fluxes, t] = self.enstrophyFluxesOverTime(timeIndices=options.timeIndices);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(1,1,TileSpacing="compact");

            for iForce = 1:length(forcing_fluxes)
                plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale,t)), hold on
            end
            legend(forcing_fluxes.fancyName)

            xlabel("time (" + self.tscale_units + ")")
            ylabel("enstrophy flux (" + self.z_flux_scale_units + ")")
            xlim([min(t) max(t)]/self.tscale);

        end

        function fig = plotEnstrophyInertialFluxOverTime(self,options)
            % Plot forcing flux for each reservoir over time
            %
            % Plots the energy flux into each reservoir from external forcing as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: fig = plotForcingFluxOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.visible: figure visibility (default: "on")
            % - Parameter options.filter: function handle to filter fluxes (default: @(v) v)
            % - Returns fig: handle to the generated figure
            arguments
                self WVDiagnostics
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
                options.timeIndices = Inf;
                options.visible = "on"
                options.filter = @(v,t) v;
            end
            [forcing_fluxes, t] = self.enstrophyInertialFluxesOverTime(timeIndices=options.timeIndices,triadComponents=options.triadComponents);

            fig = figure(Visible=options.visible);
            tl = tiledlayout(1,1,TileSpacing="compact");

            for iForce = 1:length(forcing_fluxes)
                plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale,t)), hold on
            end
            legend(forcing_fluxes.fancyName)

            xlabel("time (" + self.tscale_units + ")")
            ylabel("enstrophy flux (" + self.z_flux_scale_units + ")")
            xlim([min(t) max(t)]/self.tscale);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (spatial temporal average, from diagnostics file)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotSourcesSinksReservoirsDiagram(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Data (over time)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [reservoirs, t] = energyOverTime(self,options)
            % Compute energy for each reservoir over time
            %
            % Returns the energy in each specified reservoir as a function of time.
            %
            % - Topic: Figures (over time)
            % - Declaration: [reservoirs, t] = energyOverTime(self,options)
            % - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
            % - Parameter options.shouldIncludeExactTotalEnergy: include exact total energy (default: true)
            % - Parameter options.timeIndices: indices for time selection (default: Inf)
            % - Returns reservoirs: struct array with energy for each reservoir
            % - Returns t: time vector
            arguments
                self WVDiagnostics
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
                options.shouldIncludeExactTotalEnergy = false
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

        function [energy, t] = exactEnergyOverTime(self, options)
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter = @(v) v;
            else
                filter = @(v) v(options.timeIndices);
            end
            [ke,ape] =self.diagfile.readVariables('ke','ape');
            energy = filter(ke+ape);
            t = filter(self.diagfile.readVariables('t'));
        end

        function [enstrophy, t] = enstrophyQGPVOverTime(self, options)
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter = @(v) v;
            else
                filter = @(v) v(options.timeIndices);
            end
            enstrophy = filter(self.diagfile.readVariables('enstrophy_quadratic'));
            t = filter(self.diagfile.readVariables('t'));
        end

        function [enstrophy, t] = enstrophyAPVOverTime(self, options)
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter = @(v) v;
            else
                filter = @(v) v(options.timeIndices);
            end
            enstrophy = filter(self.diagfile.readVariables('enstrophy_apv'));
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

        function forcing_fluxes = exactForcingFluxesSpatialTemporalAverage(self,options)
            % Compute spatial-temporal average of the exact forcing fluxes
            %
            % Returns the spatial-temporal average of the exact energy fluxes from external forcing
            %
            % - Topic: Flux averages, scalar
            % - Declaration: forcing_fluxes = exactForcingFluxesSpatialTemporalAverage(self,options)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns forcing_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end

            forcing_fluxes = self.exactForcingFluxesOverTime(timeIndices=options.timeIndices);
            for iForce = 1:length(forcing_fluxes)
                forcing_fluxes(iForce).te = mean(forcing_fluxes(iForce).te);
            end
        end

        function enstrophy_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
            % Compute spatial-temporal average of the exact enstrophy fluxes
            %
            % Returns the spatial-temporal average of the exact enstrophy fluxes from external forcing
            %
            % - Topic: Flux averages, scalar
            % - Declaration: forcing_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns forcing_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end

            enstrophy_fluxes = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
            for iForce = 1:length(enstrophy_fluxes)
                enstrophy_fluxes(iForce).Z0 = mean(enstrophy_fluxes(iForce).Z0);
            end
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

        function [forcing_fluxes, t] = exactForcingFluxesOverTime(self,options)
            % Compute exact forcing fluxes over time
            %
            % Returns the exact energy fluxes from external forcing for each time step.
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: forcing_fluxes = exactForcingFluxesOverTime(self)
            % - Returns forcing_fluxes: struct array with exact fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter_space = @(v) v;
            else
                filter_space = @(v) v(options.timeIndices);
            end
            forcingNames = self.wvt.forcingNames;
            forcing_fluxes(length(forcingNames)) = struct("name","placeholder");

            for iForce=1:length(forcingNames)
                name = replace(forcingNames(iForce),"-","_");
                name = replace(name," ","_");
                forcing_fluxes(iForce).name = name;
                forcing_fluxes(iForce).fancyName = forcingNames(iForce);
                forcing_fluxes(iForce).te = filter_space(self.diagfile.readVariables("E_" + name));
            end

            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
        end

        function [forcing_fluxes, t] = exactEnstrophyFluxesOverTime(self,options)
            % Compute exact enstrophy fluxes over time
            %
            % Returns the exact enstrophy fluxes from external forcing for each time step.
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: forcing_fluxes = exactEnstrophyFluxesOverTime(self)
            % - Returns forcing_fluxes: struct array with exact fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter_space = @(v) v;
            else
                filter_space = @(v) v(options.timeIndices);
            end
            forcingNames = self.wvt.forcingNames;
            forcing_fluxes(length(forcingNames)) = struct("name","placeholder");

            for iForce=1:length(forcingNames)
                name = replace(forcingNames(iForce),"-","_");
                name = replace(name," ","_");
                forcing_fluxes(iForce).name = name;
                forcing_fluxes(iForce).fancyName = forcingNames(iForce);
                forcing_fluxes(iForce).Z0 = filter_space(self.diagfile.readVariables("Z_" + name));
            end

            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
        end

        function [forcing_fluxes,t] = forcingFluxesOverTime(self,options)
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
            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
        end

        function [inertial_fluxes,t] = inertialFluxesOverTime(self,options)
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
                options.timeIndices = Inf;
            end
            if isinf(options.timeIndices)
                filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
            else
                filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
            end
            inertial_fluxes = self.filterFluxesForReservoir(self.inertialFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
        end

        function [enstrophy_fluxes,t] = enstrophyFluxesOverTime(self,options)
            % Compute enstrophy fluxes over time
            %
            % Returns the enstrophy fluxes from external forcing 
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: enstrophy_fluxes = enstrophyFluxesOverTime(self,options)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns enstrophy_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end

            if isinf(options.timeIndices)
                filter = @(v) reshape( sum(sum(v,1),2), [], 1);
            else
                filter = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
            end
            enstrophy_fluxes = self.enstrophyFluxes;
            for iForce=1:length(enstrophy_fluxes)
                    enstrophy_fluxes(iForce).Z0 = filter(enstrophy_fluxes(iForce).Z0);
            end
            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
        end

        function [enstrophy_fluxes,t] = enstrophyInertialFluxesOverTime(self,options)
            % Compute enstrophy inertial (aka, triad) fluxes over time
            %
            % Returns the enstrophy fluxes from external forcing
            %
            % - Topic: Fluxes over time, [t 1]
            % - Declaration: enstrophy_fluxes = enstrophyFluxesOverTime(self,options)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns enstrophy_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
                options.timeIndices = Inf;
            end

            if isinf(options.timeIndices)
                filter = @(v) reshape( sum(sum(v,1),2), [], 1);
            else
                filter = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
            end
            enstrophy_fluxes = self.enstrophyInertialFluxes(triadComponents=options.triadComponents);
            for iForce=1:length(enstrophy_fluxes)
                enstrophy_fluxes(iForce).Z0 = filter(enstrophy_fluxes(iForce).Z0);
            end
            t = self.t_diag;
            if ~isinf(options.timeIndices)
                t = t(options.timeIndices);
            end
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

        function enstrophy_fluxes = enstrophyFluxesTemporalAverage(self,options)
            % Compute temporally averaged enstrophy fluxes
            %
            % Returns the temporally averaged enstrophy fluxes from external forcing for each reservoir.
            %
            % - Topic: Fluxes in space, [j kRadial]
            % - Declaration: enstrophy_fluxes = enstrophyFluxesTemporalAverage(self,options)
            % - Parameter options.timeIndices: indices for time averaging (default: Inf)
            % - Returns enstrophy_fluxes: struct array with averaged fluxes
            arguments
                self WVDiagnostics
                options.timeIndices = Inf;
            end

            if isinf(options.timeIndices)
                filter_space = @(v) mean(v,3);
            else
                filter_space = @(v) mean(v(:,:,options.timeIndices),3);
            end
            enstrophy_fluxes = self.enstrophyFluxes;
            enstrophy_fluxes.Z0 = filter_space(enstrophy_fluxes.Z0);
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
                fluxes = EnergyReservoir.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
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
                for m=1:length(primaryFlowComponents_t)
                    name = primaryFlowComponents_t(i).abbreviatedName + "_" + primaryFlowComponents_t(m).abbreviatedName;

                    % these are temporary variavbles for use within this loop only
                    Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
                    Ejk.Em = self.diagfile.readVariables("Em_" + name);
                    Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
                    Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

                    % total fluxes
                    inertial_fluxes(counter).name = name;
                    inertial_fluxes(counter).fancyName = primaryFlowComponents_t(i).shortName + "-" + primaryFlowComponents_t(m).shortName;

                    % per-reservoir fluxes
                    fluxes = EnergyReservoir.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
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
                for m=1:n
                    counter = m+(i-1)*n;
                    inertial_fluxes_consol(counter).name = options.triadComponents(i).name + "_" + options.triadComponents(m).name;
                    inertial_fluxes_consol(counter).fancyName = options.triadComponents(i).name + "{\nabla}" + options.triadComponents(m).name;
                    for iReservoir = 1:length(options.energyReservoirs)
                        arraySize = size(inertial_fluxes(1).(options.energyReservoirs(iReservoir).name));
                        inertial_fluxes_consol(counter).(options.energyReservoirs(iReservoir).name) = zeros(arraySize);
                    end
                end
            end

            counter = 1;
            for i=1:length(primaryFlowComponents_t)
                index_i = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(i)), consolidatedFlowComponents_t));
                for m=1:length(primaryFlowComponents_t)
                    index_j = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(m)), consolidatedFlowComponents_t));
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

        function enstrophy_fluxes = enstrophyFluxes(self)
            % Return the enstrophy flux from the forcing terms
            %
            % Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
            %
            % - Topic: Core function — spatial temporal
            % - Declaration: forcing_fluxes = forcingFluxes(options)
            % - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
            % - Returns forcing_fluxes: an array of structs
            arguments
                self WVDiagnostics
            end
            forcingNames = self.wvt.forcingNames;
            enstrophy_fluxes(length(forcingNames)) = struct("name","placeholder");

            for iForce=1:length(forcingNames)
                name = replace(forcingNames(iForce),"-","_");
                name = replace(name," ","_");

                enstrophy_fluxes(iForce).name = name;
                enstrophy_fluxes(iForce).fancyName = forcingNames(iForce);
                enstrophy_fluxes(iForce).Z0 = self.diagfile.readVariables("Z0_" + name);
            end
        end

        function inertial_fluxes = enstrophyInertialFluxes(self,options)
            % Return the enstrophy flux from the forcing terms
            %
            % Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
            %
            % - Topic: Core function — spatial temporal
            % - Declaration: forcing_fluxes = forcingFluxes(options)
            % - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
            % - Returns forcing_fluxes: an array of structs
            arguments
                self WVDiagnostics
                options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
            end
            primaryFlowComponents_t = self.wvt.primaryFlowComponents;
            inertial_fluxes(length(primaryFlowComponents_t)^2) = struct("name","placeholder");

            counter = 1;
            for i=1:length(primaryFlowComponents_t)
                for m=1:length(primaryFlowComponents_t)
                    name = primaryFlowComponents_t(i).abbreviatedName + "_" + primaryFlowComponents_t(m).abbreviatedName;

                    % total fluxes
                    inertial_fluxes(counter).name = name;
                    inertial_fluxes(counter).fancyName = primaryFlowComponents_t(i).shortName + "-" + primaryFlowComponents_t(m).shortName;
                    inertial_fluxes(counter).Z0 = self.diagfile.readVariables("Z0_" + name);

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
                for m=1:n
                    counter = m+(i-1)*n;
                    inertial_fluxes_consol(counter).name = options.triadComponents(i).name + "_" + options.triadComponents(m).name;
                    inertial_fluxes_consol(counter).fancyName = options.triadComponents(i).name + "{\nabla}" + options.triadComponents(m).name;
                    inertial_fluxes_consol(counter).Z0 = zeros(size(inertial_fluxes(1).Z0));
                end
            end

            counter = 1;
            for i=1:length(primaryFlowComponents_t)
                index_i = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(i)), consolidatedFlowComponents_t));
                for m=1:length(primaryFlowComponents_t)
                    index_j = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(m)), consolidatedFlowComponents_t));
                    a = inertial_fluxes(counter).Z0;
                    b = inertial_fluxes_consol(index_j+(index_i-1)*n).Z0;
                    inertial_fluxes_consol(index_j+(index_i-1)*n).Z0 = a + b;
                    counter = counter+1;
                end
            end

            inertial_fluxes = inertial_fluxes_consol;
        end

        function setLogWavelengthXAxis(self,options)
            arguments
                self WVDiagnostics
                options.num_ticks = 6
                options.roundToNearest = 5
            end
            [labels_x,ticks_x] = self.logWavelengthAxis(num_ticks=options.num_ticks,roundToNearest=options.roundToNearest);
            xscale('log')
            xticks(ticks_x)
            xticklabels(labels_x)
        end

        function [labels, ticks] = logWavelengthAxis(self,options)
            % To use this:
            % xticks(ticks_x)
            % xticklabels(labels_x)
            arguments
                self WVDiagnostics
                options.num_ticks = 6
                options.roundToNearest = 5
            end
            ticks = logspace(log10(self.kRadial(2)),log10(self.kRadial(end)),options.num_ticks);
            ticks = round(2*pi./(1e3.*ticks)/options.roundToNearest)*options.roundToNearest;
            labels = cell(length(ticks),1);
            for i=1:length(ticks)
                labels{i} = sprintf('%.0f',ticks(i));
            end
            ticks = 2*pi./(1e3*ticks);
        end

        function overlayFrequencyContours(self,options)
            arguments
                self 
                options.frequencies = [1.01 1.05 1.2 1.5 2 4 8 16]
                options.textColor = [.5,.5,.5]
                options.labelSpacing = 400
                options.lineWidth = 1
            end
            [omegaN,n] = self.wvt.transformToRadialWavenumber(abs(self.wvt.Omega),ones(size(self.wvt.Omega)));
            omegaJK = (omegaN./n)/self.wvt.f;
            set(gca,'layer','top'),
            hold on
            [C,h] = contour(2*pi./self.kRadial(2:end)/1000,self.j',(omegaJK(:,2:end)),options.frequencies,'LineWidth',options.lineWidth,'Color',options.textColor);
            clabel(C,h,options.frequencies,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
        end

        function overlayGeostrophicKineticPotentialRatioContours(self,options)
            arguments
                self 
                options.ratios = [-2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5]
                options.textColor = [.5,.5,.5]
                options.labelSpacing = 400
                options.lineWidth = 1
            end
            hke = self.wvt.transformToRadialWavenumber( self.wvt.A0_KE_factor );
            pe = self.wvt.transformToRadialWavenumber( self.wvt.A0_PE_factor );
            ratio = log10(hke./pe);
            set(gca,'layer','top'),
            hold on
            [C,h] = contour(self.kRadial(2:end),self.j(2:end)',(ratio(2:end,2:end)),options.ratios,'LineWidth',options.lineWidth,'Color',options.textColor);
            clabel(C,h,options.ratios,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
        end

        function showRossbyRadiusYAxis(self,options)
            arguments
                self 
                options.textColor = [.5,.5,.5] 
            end
            set(gca,'Layer','top','TickLength',[0.015 0.015])
            % create some nice tick labels to show deformation radius
            yticksTemp = yticks;
            ticks_y = sqrt(self.wvt.Lr2)./1000;
            labels_y = cell(length(yticksTemp),1);
            for i=1:length(yticksTemp)
                labels_y{i} = sprintf('%0.1f',ticks_y(yticksTemp(i)+1));
            end
            text(.7*min(xlim)*ones(size(yticksTemp)),yticksTemp,labels_y,'Color',options.textColor,'HorizontalAlignment','center')
            text(.7*min(xlim),1.1*max(ylim),'L_r (km)','Color',options.textColor,'HorizontalAlignment','center')
        end
    end

    methods (Static)
        cmap = cmocean(ColormapName,varargin)
        cmap = crameri(ColormapName,varargin)

        function iTimeChanged(~,eventData)
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


    end
end