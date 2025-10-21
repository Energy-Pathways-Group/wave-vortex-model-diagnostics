function fig = plotSourcesSinksForReservoirGroup(self,options)
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
    options.name string = "reservoir-damped-wave-geo"
    options.customNames = configureDictionary("string","string")
    options.fluxTolerance = 1e-2;
    options.shouldShowUnits = true;
    options.timeIndices = Inf;
    options.shouldShowReservoirEnergy = true
    options.shouldShowExactValues = true
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

[inertial, forcing, ddt, reservoir_energy] = self.fluxesForReservoirGroup(name=options.name,timeIndices=options.timeIndices,outputFormat="struct");

% get rid of the advection term
forcing(1) = [];

% add the totals for each forcing
unknownFields = setdiff(fieldnames(forcing), {'name','fancyName'});
totalForcingFlux = zeros(length(forcing),1);
for f = unknownFields'
    totalForcingFlux = totalForcingFlux + [forcing.(f{1})].';
end

% now sort into sources and sincs
iSink = 1; iSource = 1;
for iForce=1:length(forcing)
    forcing(iForce).total = totalForcingFlux(iForce);
    if totalForcingFlux(iForce) < 0
        forcing_sinks(iSink) = forcing(iForce);
        iSink = iSink + 1;
    else
        forcing_sources(iSource) = forcing(iForce);
        iSource = iSource + 1;
    end
end

reservoirs = configureDictionary("string","Box");
for iReservoir = 1:length(inertial)
    name = inertial(iReservoir).name;

    if isKey(options.customNames,name)
        fancyName = options.customNames(name);
    else
        fancyName = inertial(iReservoir).fancyName;
    end

    % reservoirs(name) = Box(fancyName,FaceColor=col{name}, FontSize=16, CornerRadius=0.10);
    reservoirs(name) = Box(fancyName, FontSize=16, CornerRadius=0.10);
    if options.shouldShowReservoirEnergy
        energy = reservoir_energy.(name)/self.escale;
        flux = ddt.(name)/self.flux_scale;

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
        if abs(forcing_fluxes(iFlux).total/self.flux_scale/2) > options.fluxTolerance
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
for i=2:length(reservoirNames)
    sourceName = inertial(i).name;
    for k=1:(i-1)
        destinationName = reservoirNames(k);
        flux = inertial(i).(destinationName)/self.flux_scale;
        magnitude = abs(flux);
        if options.shouldShowUnits
            label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
        else
            label = sprintf("%.2f",magnitude);
        end
        if magnitude > options.fluxTolerance
            if flux > 0
                inertial_arrows(end+1) = Arrow(reservoirs(sourceName),reservoirs(destinationName),Label=label,Magnitude=magnitude, LabelOffset=0.3, FontSize=14);
            else
                inertial_arrows(end+1) = Arrow(reservoirs(destinationName),reservoirs(sourceName),Label=label,Magnitude=magnitude, LabelOffset=0.3, FontSize=14);
            end
        end
    end
end

RowSublabels = strings(1,3);
if options.shouldShowExactValues
    source_total = sum([forcing_sources.te_exact]) +  sum([forcing_sources.te_exact_damp]);
    sink_total = sum([forcing_sinks.te_exact]) + sum([forcing_sinks.te_exact_damp]);

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
% fig = plotThreePointFiveRowBoxDiagram(sources, boxes(1:2), boxes(3), sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[4.5 1.5], Title=options.title, RowSublabels=RowSublabels, visible=options.visible);
fig = plotThreeRowBoxDiagram(sources, boxes, sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[4.5 1.5], Title=options.title, RowSublabels=RowSublabels, visible=options.visible);


end