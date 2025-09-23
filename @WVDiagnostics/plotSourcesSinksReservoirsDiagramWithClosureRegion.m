function fig = plotSourcesSinksReservoirsDiagramWithClosureRegion(self,options)
% Plot sources, sinks, and reservoirs diagram
%
% Generates a diagram showing energy sources, sinks, and reservoirs, including fluxes between them.
%
% - Topic: Figures (over time)
% - Declaration: fig = plotSourcesSinksReservoirsDiagram(self,options)
% - Parameter options.energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave])
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
    options.customNames = configureDictionary("string","string")
    options.fluxTolerance = 1e-2;
    options.shouldShowUnits = true;
    options.timeIndices = Inf;
    options.shouldShowReservoirEnergy = true
    options.shouldShowExactValues = true
    options.shouldSeparateClosureRegion = true
    options.title = "Energy Pathways";
    options.visible = "on"
end

col = configureDictionary("string","cell");
col{"source"} = [191 191 250]/255;
col{"ke_g"} = [205 253 254]/255;
col{"pe_g"} = [205 253 254]/255;
col{"te_gmda"} = [205 253 254]/255;
col{"te_wave"} = [205 253 197]/255;
col{"te_damp"} = [218, 160, 109]/255;
col{"sink"} = [245 194 193]/255;

[forcing_sources, forcing_sinks, inertial, ddt, reservoir_energy] = self.filterEnergyForSourcesSinksReservoirs(timeIndices=options.timeIndices);
[reservoirEnergy, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);



reservoirs = configureDictionary("string","Box");

energyReservoirs = cat(2,options.energyReservoirs,EnergyReservoir.damp);
for iReservoir = 1:length(energyReservoirs)
    name = energyReservoirs(iReservoir).name;
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

        energy = reservoir_energy.(name)/self.escale;
        if name == "te_damp"
            flux = 0;
        else
            flux = (reservoirEnergy(iReservoir).energy(end) - reservoirEnergy(iReservoir).energy(1))/(t(end)-t(1))/self.flux_scale;
        end
        
        if abs(flux) > options.fluxTolerance
            if flux > 0
                reservoirs(name).Sublabel=sprintf("%.2f %s (+%.2f %s)",energy,self.escale_units,abs(flux),self.flux_scale_units);
            else
                reservoirs(name).Sublabel=sprintf("%.2f %s (–%.2f %s)",energy,self.escale_units,abs(flux),self.flux_scale_units);
            end
        else
            reservoirs(name).Sublabel=sprintf("%.2f %s",energy,self.escale_units);
        end
    end
end

sources = Box.empty(0,0);
sinks = Box.empty(0,0);
source_arrows = Arrow.empty(0,0);
sink_arrows = Arrow.empty(0,0);

for i=1:2
    if i==1
        forcing_fluxes = forcing_sources;
    else
        forcing_fluxes = forcing_sinks;
    end
    for iFlux=1:length(forcing_fluxes)
        if isKey(options.customNames,forcing_fluxes(iFlux).name)
            fancyName = options.customNames(forcing_fluxes(iFlux).name);
        else
            fancyName = forcing_fluxes(iFlux).fancyName;
        end

        box = [];
        if abs((forcing_fluxes(iFlux).te_exact + forcing_fluxes(iFlux).te_exact_damp)/self.flux_scale/2) > options.fluxTolerance
            box = Box(fancyName,FaceColor=col{"source"}, FontSize=16);
            if i==1
                sources(end+1) = box;
            else
                sinks(end+1) = box;
            end
            reservoirNames = reservoirs.keys;
            for iRes=1:length(reservoirNames)
                name = reservoirNames(iRes);
                magnitude = abs(forcing_fluxes(iFlux).(name))/self.flux_scale;
                if options.shouldShowUnits
                    label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
                else
                    label = sprintf("%.2f",magnitude);
                end
                if abs(magnitude) > options.fluxTolerance
                    if i==1
                        source_arrows(end+1) = Arrow(sources(end),reservoirs(name),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    else
                        sink_arrows(end+1) = Arrow(reservoirs(name),sinks(end),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
                    end
                end
            end
        end

        if ~isempty(box) && options.shouldShowExactValues
            if options.shouldShowUnits
                box.Sublabel = sprintf("[%.2f %s]",(forcing_fluxes(iFlux).te_exact+forcing_fluxes(iFlux).te_exact_damp)/self.flux_scale,self.flux_scale_units);
            else
                box.Sublabel = sprintf("[%.2f]",(forcing_fluxes(iFlux).te_exact+forcing_fluxes(iFlux).te_exact_damp)/self.flux_scale);
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
flux = inertial(1).te_wave/self.flux_scale;
magnitude = abs(flux);
if options.shouldShowUnits
    label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
else
    label = sprintf("%.2f",magnitude);
end
if magnitude > options.fluxTolerance
    if flux > 0
        inertial_arrows(end+1) = Arrow(reservoirs("te_wave"),reservoirs("te_gmda"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    else
        inertial_arrows(end+1) = Arrow(reservoirs("te_gmda"),reservoirs("te_wave"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    end
end

flux = inertial(1).te_damp/self.flux_scale;
magnitude = abs(flux);
if options.shouldShowUnits
    label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
else
    label = sprintf("%.2f",magnitude);
end
if magnitude > options.fluxTolerance
    if flux > 0
        inertial_arrows(end+1) = Arrow(reservoirs("te_damp"),reservoirs("te_gmda"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    else
        inertial_arrows(end+1) = Arrow(reservoirs("te_gmda"),reservoirs("te_damp"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    end
end

flux = inertial(2).te_damp/self.flux_scale;
magnitude = abs(flux);
if options.shouldShowUnits
    label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
else
    label = sprintf("%.2f",magnitude);
end
if magnitude > options.fluxTolerance
    if flux > 0
        inertial_arrows(end+1) = Arrow(reservoirs("te_damp"),reservoirs("te_wave"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    else
        inertial_arrows(end+1) = Arrow(reservoirs("te_wave"),reservoirs("te_damp"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
    end
end


RowSublabels = strings(1,3);
if options.shouldShowExactValues
    source_total = sum([forcing_sources.te_exact]);
    sink_total = sum([forcing_sinks.te_exact]);

    if options.shouldShowUnits
        RowSublabels(1) = sprintf("[%.2f %s]",source_total/self.flux_scale,self.flux_scale_units);
        RowSublabels(3) = sprintf("[%.2f %s]",sink_total/self.flux_scale,self.flux_scale_units);
    else
        RowSublabels(1) = sprintf("[%.2f]",source_total/self.flux_scale);
        RowSublabels(3) = sprintf("[%.2f]",sink_total/self.flux_scale);
    end
    

    mean_energy = reservoir_energy.te_exact/self.escale;
    flux = ddt.te_exact/self.flux_scale;
    if abs(flux) > options.fluxTolerance
        if flux > 0
            RowSublabels(2)=sprintf("[%.2f %s (+%.2f %s)]",mean_energy,self.escale_units,abs(flux),self.flux_scale_units);
        else
            RowSublabels(2)=sprintf("[%.2f %s (–%.2f %s)]",mean_energy,self.escale_units,abs(flux),self.flux_scale_units);
        end
    else
        RowSublabels(2)=sprintf("[%.2f %s]",mean_energy,self.escale_units);
    end
end

% fig = plotThreeRowBoxDiagram(sources, reservoirs.values, sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[3.0 1.5], Title=options.title, RowSublabels=RowSublabels, visible=options.visible);
boxes = reservoirs.values;
fig = plotThreePointFiveRowBoxDiagram(sources, boxes(1:2), boxes(3), sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[4.5 1.5], Title=options.title, RowSublabels=RowSublabels, visible=options.visible);

end