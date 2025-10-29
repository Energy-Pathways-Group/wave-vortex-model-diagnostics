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
        pseudoRadialBinning (1,1) string {mustBeMember(pseudoRadialBinning,["k2+j2","adaptive"])} = "k2+j2"
    end

    properties (Access=private)
        wvt_aa_cache
    end

    properties (Dependent)
        t_diag
        t_wv
        j
        jWavenumber
        kRadial
        kPseudoRadial
        omegaAxis
        kePeAxis
        Lr2
        Lr2_pm
        omega_jk
        geo_hke_jk
        geo_pe_jk
        forcingNames
        wvt_aa
        diagnosticsHasExplicitAntialiasing
    end

    properties (SetObservable)
        iTime = Inf
    end

    methods
        createDiagnosticsFile(self,options)
        createWWGTriadDiagnostic(self,options)
        createDampedRegionDiagnostics(self,options)
        create1DMirrorFluxes(self,options)
        create2DMirrorFluxes(self,options)

        createReservoirGroup(self,options) % groupName, flowComponents
        [triadVar, forcingVar, energyVar] = variablesForReservoirGroup(self,options)
        addTriadFluxesForReservoirGroupAtTime(self,options) %triadVar, flowComponents, wvt, outputIndex
        addForcingFluxesForReservoirGroupAtTime(self,options) %forcingVar, flowComponents, wvt, outputIndex
        [transferFlux, forcingFlux, ddt, energy] = fluxesForReservoirGroup(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (instantaneous snapshot, from model output)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotFluidStateMultipanel(self,options)
        fig = plotFluidDecompositionMultipanel(self,options)

        fig = plotPotentialEnergySpectrum(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures for Energy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotEnergySpectrum(self,options)
        fig = plotEnergyOverTime(self,options)
        fig = plotEnergyFluxOverTime(self,options)
        fig = plotEnergyTriadFluxOverTime(self,options)
        fig = plotEnergyFluxTemporalAverage(self,options)

        fig = plotEnergyFluxes1D(self,options)

        fig = plotSourcesSinksReservoirsDiagram(self,options)
        fig = plotSourcesSinksReservoirsDiagramWithClosureRegion(self,options)
        summarizeSourcesSinksReservoirs(self,options)
        [sources, sinks, inertial_tx, inertial_cascade, ddt, energy] = filterEnergyForSourcesSinksReservoirs(self,options)

        showDampingFluxVsPseudolength(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures for Potential Enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotEnstrophySpectrum(self,options)
        fig = plotEnstrophyOverTime(self,options)
        fig = plotEnstrophyFluxOverTime(self,options)
        fig = plotEnstrophyTriadFluxOverTime(self,options)
        fig = plotEnstrophyFluxTemporalAverage(self,options)

        tableString = createEnstrophyFluxSummaryTable(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Figures (ancillary data, from model output)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotMooringRotarySpectrum(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Data (over time)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [energy, t] = exactEnergyOverTime(self, options)
        [enstrophy, t] = exactEnstrophyOverTime(self, options)

        [reservoirs, t] = quadraticEnergyOverTime(self,options)
        [enstrophy, t] = quadraticEnstrophyOverTime(self, options)

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

        [E0_ggw,Epm_wwg] = quadraticEnergyMirrorTriadsUndamped(self,options);

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
        % Fluxes in space, [sparseKRadialAxis 1]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [M_wwg, M_ggw, kp] = quadraticEnergyMirrorTriadFluxes1D(self,options)
        [M_wwg, omegaAxis] = quadraticEnergyMirrorTriadFluxes1D_omega(self,options)
        [M_ggw, kePeAxis] = quadraticEnergyMirrorTriadFluxes1D_kepe(self,options)

        [inertial_fluxes_g, inertial_fluxes_w, kp] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(self,options)
        [inertial_fluxes_w, omegaAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega(self,options)
        [inertial_fluxes_g, kePeAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe(self,options)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [M_wwg, F_wwg, ks, js] = quadraticEnergyMirrorTriadFluxes2D(self,options)

        [inertial_fluxes_g, inertial_fluxes_w, ks, js] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)

        wwg = quadraticEnergyTriadFluxWWGWave(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transforms to alternative axes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [varargout] = transformToPseudoRadialWavenumber(self,energyReservoir,varargin);
        [varargout] = transformToPseudoRadialWavenumberA0(self,varargin);
        [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)  
        [varargout] = transformToOmegaAxis(self,varargin) 
        [varargout] = transformToKePeAxis(self,varargin) 

        [kp,bins_0,bins_pm] = sparsePseudoRadialAxis(self)
        [omegaAxis,bins_omega] = sparseOmegaAxis(self)
        [kePeAxis,bins_kepe] = sparseKePeAxis(self)
        k = sparseKRadialAxis(self)
        j = sparseJWavenumberAxis(self)
        [S_0, S_pm, mask_0, mask_pm] = sparseJKAxisBinMatrices(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Spectra (also implemented in the WVT)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S_f = spectrumWithFgTransform(self,f)
        S_f = spectrumWithGgTransform(self,f)
        S_f = crossSpectrumWithFgTransform(self,phi,gamma)
        S_f = crossSpectrumWithGgTransform(self,phi,gamma)

        function fig = plotPoissonFlowOverPcolor(wvd,options)
            arguments
                wvd WVDiagnostics
                options.visible = "on"
                options.vectorFlux
                options.flux
                options.overSaturationFactor = 2;
            end

            color_axis_limits = max(abs(options.flux(:)))*[-1 1]/options.overSaturationFactor;

            % We will pretend the "0" wavenumber is actually evenly spaced
            % from the nearest two wavenumbers
            kPseudoLocation = wvd.kRadial;
            kPseudoLocation(1) = exp(-log(kPseudoLocation(3)) + 2*log(kPseudoLocation(2)));
            jPseudoLocation = wvd.jWavenumber;
            jPseudoLocation(1) = exp(-log(jPseudoLocation(3)) + 2*log(jPseudoLocation(2)));

            % For interpolation to work correctly we need to repeat the
            % first entry, but properly back at zero
            kPadded = cat(1,0,kPseudoLocation);
            jPadded= cat(1,0,jPseudoLocation);
            [KPadded,JPadded] = ndgrid(kPadded,jPadded);
            fluxPadded = cat(1,options.flux(1,:),options.flux);
            fluxPadded = cat(2,fluxPadded(:,1),fluxPadded);

            % We will use this axis to display. Widen the box so that it is
            % the same size as its neighbors.
            kMin = exp(-1.5*log(wvd.kRadial(3)) + 2.5*log(wvd.kRadial(2)));
            jMin = exp(-1.5*log(wvd.jWavenumber(3)) + 2.5*log(wvd.jWavenumber(2)));

            N = 500;
            kLinLog = linspace(log10(kMin),max(log10(wvd.kRadial(:))),N);
            jLinLog = linspace(log10(jMin),max(log10(wvd.jWavenumber(:))),N/2);
            [KLinLog,JLinLog] = ndgrid(kLinLog,jLinLog);
            fluxLinLog = interpn(KPadded,JPadded,(fluxPadded.'),10.^KLinLog,10.^JLinLog,"nearest");

            fig = figure(Units='points',Visible = options.visible);
            set(gcf,'PaperPositionMode','auto')

            pcolor(KLinLog,JLinLog,fluxLinLog); shading flat;
            colormap(WVDiagnostics.crameri('-bam'))
            clim(color_axis_limits)
            colorbar("eastoutside")

            [X,Y,U,V] = wvd.PoissonFlowFromFlux(options.vectorFlux.');
            [logX,logY,Uprime,Vprime] = wvd.RescalePoissonFlowFluxForLogSpace(X,Y,U,V,shouldOnlyRescaleDirection=false);
            % we need two adjustments. First, we need to move the first row and column
            % half an increment
            logX(1,:) = log10(kPseudoLocation(1));
            logY(:,1) = log10(jPseudoLocation(1));
            scale = 1.5;
            hold on
            % quiver(logX,logY,scale*Uprime,scale*Vprime,'off',Color=0*[1 1 1],AutoScale="off",LineWidth=2)
            quiver(logX,logY,Uprime,Vprime,Color=0*[1 1 1],AutoScale="off",LineWidth=2)
        end

fig = plotPoissonFlowOverContours(wvd,options)

[logX,logY,Uprime,Vprime] = RescalePoissonFlowFluxForLogSpace(wvd,X,Y,U,V,options)

        function [X,Y,U,V] = PoissonFlowFromFlux(wvd, flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
            % For the DCT2/DST2 we use a half-shift grid
            x = wvd.kRadial + 0*(wvd.kRadial(2)-wvd.kRadial(1))/2;
            y = wvd.jWavenumber + 0*(wvd.jWavenumber(2)-wvd.jWavenumber(1))/2;
            [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFluxWithAxes(x,y,flux);

            U = U/wvd.kRadial(2);
            V = V/wvd.jWavenumber(2);
        end

        function [X,Y,U,V] = PoissonFlowFromFluxType1(wvd, flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
            % For the DCT2/DST2 we use a half-shift grid
            x = wvd.kRadial;
            y = wvd.jWavenumber;

            N = length(x);
            M = length(y);
            dk = 1/(2*(N-1)*(x(2)-x(1)));
            dl = 1/(2*(M-1)*(y(2)-y(1)));

            k=dk*(0:(N-1)).';
            l=dl*(0:(M-1)).';

            DCTx = WVDiagnostics.DCT1(N);
            DCTy = WVDiagnostics.DCT1(M);

            flux_ky = DCTx*flux;
            flux_kl = shiftdim(DCTy*shiftdim(flux_ky,1),1);

            [K,L] = ndgrid(k,l);

            D = -((2*pi*K).^2 + (2*pi*L).^2);
            D(1,1) = Inf;

            UFactor = 2*pi*K./D;
            VFactor = 2*pi*L./D;

            iDCTx = WVDiagnostics.iDCT1(N);
            iDSTy = WVDiagnostics.iDST1(M);
            iDSTy = cat(2,zeros(M,1),iDSTy,zeros(M,1));

            V_xl = iDCTx*(VFactor.*flux_kl);
            V = shiftdim(iDSTy*shiftdim(V_xl,1),1);

            iDCTy = WVDiagnostics.iDCT1(M);
            iDSTx = WVDiagnostics.iDST1(N);
            iDSTx = cat(2,zeros(N,1),iDSTx,zeros(N,1));

            U_xl = iDSTx*(UFactor.*flux_kl);
            U = shiftdim(iDCTy*shiftdim(U_xl,1),1);

            [X,Y] = ndgrid(x,y);
        end

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
            omegaJK = self.omega_jk;
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
            hke = self.geo_hke_jk;
            pe = self.geo_pe_jk;
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
            hke = self.geo_hke_jk;
            pe = self.geo_pe_jk;
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dependent property implementations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function wvt_aa = get.wvt_aa(self)
            if isempty(self.wvt_aa_cache) || ~isvalid(self.wvt_aa_cache)
                if isa(self.wvt,"WVTransformBoussinesq")
                    if ~exist(self.wvaapath,"file")
                        fprintf("No existing Boussinesq transform with explicity antialiasing found. Creating a new transform.\n");
                        wvt_ = self.wvt.waveVortexTransformWithExplicitAntialiasing();
                        ncfile_ = wvt_.writeToFile(self.wvaapath);
                        ncfile_.close;
                    end
                    self.wvt_aa_cache = WVTransform.waveVortexTransformFromFile(self.wvaapath,iTime=self.iTime,shouldReadOnly=true);
                else
                    self.wvt_aa_cache = self.wvt.waveVortexTransformWithExplicitAntialiasing();
                end
                self.wvt_aa_cache.addOperation(EtaTrueOperation());
                self.wvt_aa_cache.addOperation(APEOperation(self.wvt_aa_cache));
                self.wvt_aa_cache.addOperation(APVOperation());
                self.wvt_aa_cache.addOperation(SpatialForcingOperation(self.wvt_aa_cache));
            end
            wvt_aa = self.wvt_aa_cache;
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

        function t = get.Lr2_pm(self)
            % Wave deformation radius squared
            %
            % Reads the 'Lr2_pm' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.Lr2_pm(self)
            % - Returns t: Lr2_pm
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('Lr2_pm');
            else
                t = self.wvt.g*self.wvt.h_pm/self.wvt.f/self.wvt.f;
            end
        end

        function t = get.omega_jk(self)
            % Wave frequency omega in jk space
            %
            % Reads the 'omega_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.omega_jk(self)
            % - Returns t: omega_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('omega_jk')
                    t = self.diagfile.readVariables('omega_jk');
                else
                    t = self.diagfile.readVariables('omega_axis');
                end
            else
                [omegaN,n] = self.wvt.transformToRadialWavenumber(abs(self.wvt.Omega),ones(size(self.wvt.Omega)));
                t = (omegaN./n);
            end
        end

        function t = get.geo_hke_jk(self)
            % Geostrophic kinetic energy in jk space
            %
            % Reads the 'geo_hke_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.geo_hke_jk(self)
            % - Returns t: geo_hke_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('geo_hke_jk')
                    t = self.diagfile.readVariables('geo_hke_jk');
                else
                    t = self.diagfile.readVariables('geo_hke_axis');
                end
            else
                t = self.wvt.transformToRadialWavenumber(self.wvt.A0_KE_factor);
            end
        end

        function t = get.geo_pe_jk(self)
            % Geostrophic potential energy in jk space
            %
            % Reads the 'geo_pe_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.geo_pe_jk(self)
            % - Returns t: geo_pe_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('geo_pe_jk')
                    t = self.diagfile.readVariables('geo_pe_jk');
                else
                    t = self.diagfile.readVariables('geo_pe_axis');
                end
            else
                t = self.wvt.transformToRadialWavenumber(self.wvt.A0_PE_factor);
            end
        end

        function jWavenumber = get.jWavenumber(self)
            jWavenumber = 1./sqrt(self.Lr2);
            jWavenumber(1) = 0; % barotropic mode is a mean?
        end

        function kPseudoRadial = get.kPseudoRadial(self)
            if self.pseudoRadialBinning == "k2+j2"
                [kj,kr] = ndgrid(self.jWavenumber,self.kRadial);
                Kh = sqrt(kj.^2 + kr.^2);
                allKs = unique(reshape(abs(Kh),[],1),'sorted');
                deltaK = max(diff(allKs));

                kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);
            elseif self.pseudoRadialBinning == "adaptive"
                % n = find(self.kRadial < self.jWavenumber(2),1,"last");
                n = find(self.kRadial > self.jWavenumber(2),1,"first");
                kAxis_ = cat(1,self.kRadial(1:n),self.jWavenumber(3:end));
            end

            kPseudoRadial = reshape(kAxis_,[],1);

            % kPseudoRadial = reshape(kAxis_,[],1);
            % 
            % n = find(kPseudoRadial < self.jWavenumber(2),1,"last");
            % kPseudoRadial = cat(1,kPseudoRadial(1:n),self.jWavenumber(2:end));

            % n = find(self.kRadial < self.jWavenumber(2),1,"last");
            % kPseudoRadial = cat(1,self.kRadial(1:n),self.jWavenumber(2:end));

        end

        function omega = get.omegaAxis(self)
            omega=self.wvt.Omega(2:end,:);
            omegaj1=omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));
            omega=min(omega(:)):dOmega:max(omega(:));
        end

        function kepe = get.kePeAxis(self)
            % fraction = self.geo_hke_jk./(self.geo_hke_jk+self.geo_pe_jk);
            a = cat(1,0,10.^linspace(log10(0.01),log10(0.5),10).');
            kepe = cat(1,a,1-flip(a(1:end-1)));
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

        function forcingNames = get.forcingNames(self)
            % Get names of the forcings
            %
            % Reads the 'forcingNames' from the wave-vortex file and adds
            % in antialiasing if appropriate.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.forcingNames(self)
            % - Returns forcingNames: strings
            forcingNames = self.wvt.forcingNames;
            if self.diagfile.hasVariableWithName("Z0_antialias_filter")
                forcingNames(end+1) = "antialias filter";
            end
        end   

        function bool = get.diagnosticsHasExplicitAntialiasing(self)
            bool = false;
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName("Z0_antialias_filter")
                    bool = true;
                end
            end
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


            [self.wvt, self.wvfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=Inf,shouldReadOnly=true);
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
    end

    methods (Static)
        cmap = cmocean(ColormapName,varargin)
        cmap = crameri(ColormapName,varargin)

        [x, t] = CosineTransformBack( f, xbar, varargin )
        [xbar, f] = CosineTransformForward( t, x, varargin )
        [x, t] = SineTransformBack( f, xbar, varargin )

        matrix = DCT1(N)
        matrix = iDCT1(N)
        matrix = DST1(N)
        matrix = iDST1(N)

        matrix = DCT2(N)
        matrix = iDCT2(N)
        matrix = DST2(N)
        matrix = iDST2(N)

        Epm = geostrophicGeostrophicWaveEnergy(wvt,mask)
        E0 = waveWaveGeostrophicEnergy(wvt,mask)
        E0 = waveWaveGeostrophicEnergyForMode(wvt,maskKU,maskKUx,Nj)

        function cmap = symmetricTintMap(c,options)
            %SYMMETRICTINTMAP Colormap going from color -> tinted white -> color
            %
            %   c   : 1x3 RGB color in [0 1]
            %   N   : total number of levels (even recommended)
            %   tintStrength : fraction of c mixed into white at the center
            %
            % Example:
            %   c = [0.2 0.6 0.8];
            %   colormap(symmetricTintMap(c,256,0.05)); colorbar
            arguments
                c 
                options.N = 256
                options.tintStrength = 0.05;
            end

            N = options.N;
            tintStrength = options.tintStrength;

            % midpoint color: mostly white with a hint of c
            c_mid = (1 - tintStrength)*[1 1 1] + tintStrength*c;

            % split N into two halves
            n1 = floor(N/2);
            n2 = N - n1;

            % first half: c -> c_mid
            half1 = [linspace(c(1), c_mid(1), n1)' ...
                linspace(c(2), c_mid(2), n1)' ...
                linspace(c(3), c_mid(3), n1)'];

            % second half: c_mid -> c
            half2 = [linspace(c_mid(1), c(1), n2)' ...
                linspace(c_mid(2), c(2), n2)' ...
                linspace(c_mid(3), c(3), n2)'];

            cmap = [half1; half2];
        end

        function iTimeChanged(~,eventData)
            self = eventData.AffectedObject;
            self.wvt.initFromNetCDFFile(self.wvfile,iTime=self.iTime);
            if ~isempty(self.wvt_aa_cache)
                self.wvt_aa_cache.t = self.wvt.t;
                [self.wvt_aa_cache.A0,self.wvt_aa_cache.Ap,self.wvt_aa_cache.Am] = self.wvt.spectralVariableWithResolution(self.wvt_aa_cache,self.wvt.A0,self.wvt.Ap,self.wvt.Am);
            end
        end

        function [X,Y,U,V] = PoissonFlowFromFluxWithAxes(x, y, flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
            % For the DCT2/DST2 we use a half-shift grid
            N = length(x);
            M = length(y);
            dk = 1/(2*N*(x(2)-x(1)));
            dl = 1/(2*M*(y(2)-y(1)));

            k=dk*(0:(N-1)).';
            l=dl*(0:(M-1)).';

            DCTx = WVDiagnostics.DCT2(N);
            DCTy = WVDiagnostics.DCT2(M);

            flux_ky = DCTx*flux;
            flux_kl = shiftdim(DCTy*shiftdim(flux_ky,1),1);

            [K,L] = ndgrid(k,l);

            D = -((2*pi*K).^2 + (2*pi*L).^2);
            D(1,1) = Inf;

            UFactor = 2*pi*K./D;
            VFactor = 2*pi*L./D;

            iDCTx = WVDiagnostics.iDCT2(N);
            iDSTy = WVDiagnostics.iDST2(M);
            iDSTy = circshift(iDSTy,1,2);
            iDSTy(:,1) = 0;

            V_xl = iDCTx*(VFactor.*flux_kl);
            V = shiftdim(iDSTy*shiftdim(V_xl,1),1);

            iDCTy = WVDiagnostics.iDCT2(M);
            iDSTx = WVDiagnostics.iDST2(N);
            iDSTx = circshift(iDSTx,1,2);
            iDSTx(:,1) = 0;

            U_xl = iDSTx*(UFactor.*flux_kl);
            U = shiftdim(iDCTy*shiftdim(U_xl,1),1);

            [X,Y] = ndgrid(x,y);
        end

        function [X,Y,U,V] = PoissonFlowFromFluxDCTI(x,y,flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])

            [X,Y] = ndgrid(x,y);
            [flux_bar, f_alpha] = WVDiagnostics.CosineTransformForward( x, flux, 1 );
            [flux_bar2, f_beta] = WVDiagnostics.CosineTransformForward( y, flux_bar, 2 );
            [ALPHA,BETA] = ndgrid(f_alpha,f_beta);
            D = -((2*pi*ALPHA).^2 + (2*pi*BETA).^2);
            D(1,1) = Inf;
            UFactor = 2*pi*ALPHA./D;
            VFactor = 2*pi*BETA./D;
            tmp = WVDiagnostics.CosineTransformBack(f_beta,UFactor.*flux_bar2,2);
            U = WVDiagnostics.SineTransformBack(f_alpha(2:end-1,:),tmp(2:end-1,:),1);
            V = WVDiagnostics.CosineTransformBack(f_alpha,WVDiagnostics.SineTransformBack(f_beta(2:end-1),VFactor(:,2:end-1).*flux_bar2(:,2:end-1),2),1);
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
                case "te_damp"
                    fancyName = "closure region";
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