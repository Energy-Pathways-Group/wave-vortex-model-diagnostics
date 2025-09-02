classdef WVDiagnostics < handle
    %WVDiagnostics Produces diagnostics and figures from WVModel output
    %   This is a collection of diagnostic tools for analyzing model output
    properties
        wvpath
        diagpath
        wvaapath
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
        jWavenumber
        kRadial
        kPseudoRadial
        omegaAxis
        Lr2
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
                    self.diagpath= fullfile(pwd,strcat(fname,"-diagnostics.nc"));
                end
            else
                self.diagpath = options.diagnosticsFilePath;
            end
            if exist(self.diagpath,"file")
                self.diagfile = NetCDFFile(self.diagpath);
            else
                warning("No diagnostics file found. Some functionality will not be available.")
            end

            if ~isempty(fpath)
                self.wvaapath = fullfile(fpath,strcat(fname,"-wvt-aa.nc"));
            else
                self.wvaapath= fullfile(pwd,strcat(fname,"-wvt-aa.nc"));
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
            % t = self.diagfile.readVariables('j');
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('j');
            else
                t = self.wvt.j;
            end
        end

        function t = get.kRadial(self)
            % Get kRadial indices
            %
            % Reads the 'kRadial' variable from the diagnostics file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.kRadial(self)
            % - Returns t: kRadial indices from diagnostics file
            % t = self.diagfile.readVariables('kRadial');
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('kRadial');
            else
                t = self.wvt.kRadial;
            end
        end

        function t = get.Lr2(self)
            % Hydrostatic deformation radius squared
            %
            % Reads the 'Lr2' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.Lr2(self)
            % - Returns t: kRadial indices from diagnostics file
            % t = self.diagfile.readVariables('kRadial');
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('Lr2');
            else
                t = self.wvt.Lr2;
            end
        end

        function jWavenumber = get.jWavenumber(self)
            jWavenumber = 1./sqrt(self.Lr2);
            jWavenumber(1) = 0; % barotropic mode is a mean?
        end

        function kPseudoRadial = get.kPseudoRadial(self)
            [kj,kr] = ndgrid(self.jWavenumber,self.kRadial);
            Kh = sqrt(kj.^2 + kr.^2);
            allKs = unique(reshape(abs(Kh),[],1),'sorted');
            deltaK = max(diff(allKs));
            kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);
            % This choices of axis spacing ensures that there will be no
            % gaps in the resulting spectrum.
            kPseudoRadial = reshape(kAxis_,[],1);
        end

        function omega = get.omegaAxis(self)
            omega=self.wvt.Omega(2:end,:);
            omegaj1=omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));
            omega=min(omega(:)):dOmega:max(omega(:));
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
        [varargout] = transformToPseudoRadialWavenumber(self,energyReservoir,varargin);
        [varargout] = transformToPseudoRadialWavenumberA0(self,varargin);
        [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)  
        [varargout] = transformToOmegaAxis(self,varargin)   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S_f = spectrumWithFgTransform(self,f)
        S_f = spectrumWithGgTransform(self,f)
        S_f = crossSpectrumWithFgTransform(self,phi,gamma)
        S_f = crossSpectrumWithGgTransform(self,phi,gamma)

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
        
        fig = plotEnstrophyOverTime(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (over time, from diagnostics file)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotEnergyOverTime(self,options)

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
            [forcing_fluxes, t] = self.quadraticEnergyFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

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
            [inertial_fluxes,t] = self.quadraticEnergyTriadFluxesOverTime(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

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
                options.shouldShowNonlinearAdvection = true
                options.shouldShowTotal = true
                options.shouldShowDtEnstrophy = true
            end
            [forcing_fluxes, t] = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
            if ~options.shouldShowNonlinearAdvection
                forcing_fluxes(1) = [];
            end

            fig = figure(Visible=options.visible);
            tl = tiledlayout(1,1,TileSpacing="compact");
            total = zeros(size(forcing_fluxes(1).Z0));
            for iForce = 1:length(forcing_fluxes)
                plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale)), hold on
                total = total + forcing_fluxes(iForce).Z0;
            end
            if options.shouldShowTotal
                plot(t/self.tscale,options.filter(total/self.z_flux_scale),Color=0*[1 1 1],LineWidth=2), hold on
                legendValues = cat(1,forcing_fluxes.fancyName,"total");
            else
                legendValues = forcing_fluxes.fancyName;
            end

            if options.shouldShowDtEnstrophy
                Z_apv = self.enstrophyAPVOverTime(timeIndices=options.timeIndices);
                t2 = t(2:end) - (t(2)-t(1))/2;
                dZdt = diff(Z_apv)./diff(t);
                plot(t2/self.tscale,options.filter(dZdt/self.z_flux_scale),Color=0*[1 1 1],LineWidth=2,LineStyle="--"), hold on
                legendValues = cat(1,legendValues,"$\frac{d Z}{dt}$");
            end

            legend(legendValues);

            xlabel("time (" + self.tscale_units + ")")
            ylabel("flux (" + self.z_flux_scale_units + ")")
            xlim([min(t) max(t)]/self.tscale);

            mean(total/self.z_flux_scale)
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
                options.shouldShowNonlinearAdvection = false
            end
            [forcing_fluxes, t] = self.quadraticEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
            if ~options.shouldShowNonlinearAdvection
                forcing_fluxes(1) = [];
            else
                forcing_fluxes(1).Z0 = forcing_fluxes(1).Z0;
            end

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
            [forcing_fluxes, t] = self.quadraticEnstrophyTriadFluxesOverTime(timeIndices=options.timeIndices,triadComponents=options.triadComponents);

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

        fig = plotEnstrophyFluxTemporalAverage(self,options)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (spatial temporal average, from diagnostics file)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotSourcesSinksReservoirsDiagram(self,options)
        
        tableString = createEnstrophyFluxSummaryTable(self,options)

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
                options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total] 
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
        % Flux averages, scalar [1 1]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        energy_fluxes = exactEnergyFluxesSpatialTemporalAverage(self,options)
        enstrophy_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)

        forcing_fluxes = quadraticEnergyFluxesSpatialTemporalAverage(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxesSpatialTemporalAverage(self,options)

        enstrophy_fluxes = quadraticEnstrophyFluxesSpatialTemporalAverage(self,options)
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes, [j kRadial t]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        energy_fluxes = exactEnergyFluxes(self)
        enstrophy_fluxes = exactEnstrophyFluxes(self)

        forcing_fluxes = quadraticEnergyFluxes(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxes(self,options)

        enstrophy_fluxes = quadraticEnstrophyFluxes(self)
        inertial_fluxes = quadraticEnstrophyTriadFluxes(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes over time, [t 1]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [energy_fluxes, t] = exactEnergyFluxesOverTime(self,options)
        [enstrophy_fluxes, t] = exactEnstrophyFluxesOverTime(self,options)

        [energy_fluxes,t] = quadraticEnergyFluxesOverTime(self,options)
        [energy_fluxes,t] = quadraticEnergyTriadFluxesOverTime(self,options)

        [enstrophy_fluxes,t] = quadraticEnstrophyFluxesOverTime(self,options)
        [enstrophy_fluxes,t] = quadraticEnstrophyTriadFluxesOverTime(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes in space, [j kRadial]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        energy_fluxes = quadraticEnergyFluxesTemporalAverage(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxesTemporalAverage(self,options)

        energy_fluxes = exactEnergyFluxesTemporalAverage(self,options)

        enstrophy_fluxes = quadraticEnstrophyFluxesTemporalAverage(self,options)
        enstrophy_fluxes = exactEnstrophyFluxesTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figure extras
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                options.labelSpacing = 600
                options.lineWidth = 1
            end
            [omegaN,n] = self.wvt.transformToRadialWavenumber(abs(self.wvt.Omega),ones(size(self.wvt.Omega)));
            omegaJK = (omegaN./n)/self.wvt.f;
            set(gca,'layer','top'),
            hold on
            % flipud() and fliplr() help trick clabel into nicer label placement. 
            % for y-axis, use j+1 so contours line up with pcolor cells.
            [C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),self.j',fliplr(omegaJK(:,2:end)),options.frequencies,'LineWidth',options.lineWidth,'Color',options.textColor);
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
            % [C,h] = contour(self.kRadial(2:end),self.j(2:end)',(ratio(2:end,2:end)),options.ratios,'LineWidth',options.lineWidth,'Color',options.textColor);
            [C,h] = contour(2*pi./self.kRadial(2:end)/1000,self.j(2:end)',(ratio(2:end,2:end)),options.ratios,'LineWidth',options.lineWidth,'Color',options.textColor);
            clabel(C,h,options.ratios,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
            [C,h] = contour(2*pi./self.kRadial(2:end)/1000,self.j(1:end)',(ratio(1:end,2:end)),[-3,3],'LineWidth',options.lineWidth,'Color',options.textColor);
            clabel(C,h,options.ratios,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
        end

        function overlayGeostrophicKineticPotentialFractionContours(self,options)
            arguments
                self 
                options.fractions = [.01,.1,.25,.75,.9,.99]
                options.textColor = [.5,.5,.5]
                options.labelSpacing = 600
                options.lineWidth = 1
            end
            hke = self.wvt.transformToRadialWavenumber( self.wvt.A0_KE_factor );
            pe = self.wvt.transformToRadialWavenumber( self.wvt.A0_PE_factor );
            fraction = hke./(hke+pe);
            set(gca,'layer','top'),
            hold on          
            % flipud() and fliplr() help trick clabel into nicer label placement. 
            % for y-axis, use j+1 so contours line up with pcolor cells.
            [C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),self.j(1:end)'+1,fliplr(fraction(1:end,2:end)),options.fractions,'LineWidth',options.lineWidth,'Color',options.textColor);
            clabel(C,h,options.fractions,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
            [C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),(self.j(1:end)')+1,fliplr(fraction(1:end,2:end)),[.5,.5],'LineWidth',2,'Color',options.textColor);
            clabel(C,h,.5,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
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

        [x, t] = CosineTransformBack( f, xbar, varargin )
        [xbar, f] = CosineTransformForward( t, x, varargin )
        [x, t] = SineTransformBack( f, xbar, varargin )

        function iTimeChanged(~,eventData)
            wvd = eventData.AffectedObject;
            wvd.wvt.initFromNetCDFFile(wvd.wvfile,iTime=wvd.iTime);
        end

        function [X,Y,U,V] = PoissonFlowFromFlux(x,y,flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])

            [X,Y] = ndgrid(x,y);
            [flux_bar, f_alpha] = CosineTransformForward( x, flux, 1 );
            [flux_bar2, f_beta] = CosineTransformForward( y, flux_bar, 2 );
            [ALPHA,BETA] = ndgrid(f_alpha,f_beta);
            D = -((2*pi*ALPHA).^2 + (2*pi*BETA).^2);
            D(1,1) = Inf;
            UFactor = 2*pi*ALPHA./D;
            VFactor = 2*pi*BETA./D;
            tmp = CosineTransformBack(f_beta,UFactor.*flux_bar2,2);
            U = SineTransformBack(f_alpha(2:end-1,:),tmp(2:end-1,:),1);
            V = CosineTransformBack(f_alpha,SineTransformBack(f_beta(2:end-1),VFactor(:,2:end-1).*flux_bar2(:,2:end-1),2),1);
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