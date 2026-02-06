---
layout: default
title: Energy flow box diagram
parent: Tutorials
mathjax: true
nav_order: 2
---

#  Customizing energy flow box diagrams

Here is a somewhat complicated example of creating an energy flow box diagram. This example shows a way to customize the spacing of the reservoir boxes, adjust the path of arrows between reservoir boxes, and adjust the placement of labels on the arrows.

```matlab
% load model diagnostics
filename = fullfile('fullpath','WaveVortexModelOutput.nc');
wvd = WVDiagnostics(filename);

% setup forcing
custom_names = configureDictionary("string","cell");
custom_names{"quadratic_bottom_friction"} = "bottom friction";
custom_names{"vertical_diffusivity"} = "diffusivity";
custom_names{"adaptive_damping"} = "damping";
custom_names{"inertial_forcing"} = "NIO";
custom_names{"M2_tidal_forcing"} = "M2 tide";
custom_names{"geostrophic_mean_flow"} = "mean flow";
custom_names{"damped_geostrophic"} = ["damped", "geostrophic"];
custom_names{"damped_wave"} = ["damped", "wave"];
customForcing = ["quadratic_bottom_friction", "adaptive_damping", "inertial_forcing", "M2_tidal_forcing","geostrophic_mean_flow"];

% setup colors
col = configureDictionary("string","cell");
col{"source"} = [191 191 250]/255;
col{"damped_geostrophic"} = [205 253 254]/255;
col{"geostrophic"} = [205 253 254]/255;
col{"wave"} = [205 253 197]/255;
col{"damped_wave"} = [205 253 197]/255;
col{"sink"} = [245 194 193]/255;

% initialize the box diagram with default spacing
order = ["geostrophic", "wave", "damped_geostrophic", "damped_wave"];
fluxTolerance = 3e-2;
[~, boxDiagram] = wvd.plotSourcesSinksForReservoirGroup(customForcing=customForcing,customNames=custom_names,customColors=col,customReservoirOrder=order,shouldShowUnits=true,fluxTolerance=fluxTolerance,timeIndices=timeIndices,title="",visible="on");

% fine tune box spacing and arrows
sourcesBoxes = boxDiagram.rows{1};
inertialBoxes = boxDiagram.rows{2};
dampedBoxes = boxDiagram.rows{3};
dampedBoxes(1).Size = [3.5 2.0];
dampedBoxes(2).Size = [3.5 2.0];
dampedBoxes(1).Position = dampedBoxes(1).Position + [4.5 0];
dampedBoxes(2).Position = dampedBoxes(2).Position + [0.5 0];
% adjust arrows and labels
sourceLabels = arrayfun(@(a) join(string(a.Source.Label), newline), boxDiagram.arrows);
targetLabels = arrayfun(@(a) join(string(a.Target.Label), newline), boxDiagram.arrows); 
% Arrow: geostrophic to bottom friction
idx = find(strcmp(sourceLabels,"geostrophic") & strcmp(targetLabels,"bottom friction"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.4;
end
% Arrow: geostrophic to wave
idx = find(strcmp(sourceLabels,"geostrophic") & strcmp(targetLabels,"wave"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.4;
end
% Arrow: geostrophic to damped wave
idx = find(strcmp(sourceLabels,"geostrophic") & strcmp(targetLabels,"damped"+newline+"wave"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.33;
end
% Arrow: NIO to damped wave
idx = find(strcmp(sourceLabels,"NIO") & strcmp(targetLabels,"damped"+newline+"wave"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.33;
end
% Arrow: wave to damped wave
idx = find(strcmp(sourceLabels,"wave") & strcmp(targetLabels,"damped"+newline+"wave"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.33;
end
% Arrow: wave to damped geostrophic
idx = find(strcmp(sourceLabels,"wave") & strcmp(targetLabels,"damped"+newline+"geostrophic"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.33;
end
% Arrow: damped geostrophic to damping
idx = find(strcmp(sourceLabels,"damped"+newline+"geostrophic") & strcmp(targetLabels,"damping"));
if ~isempty(idx)
    boxDiagram.arrows(idx).LabelOffset = 0.33;
end
% Arrow: wave to bottom friction
idx = find(strcmp(sourceLabels,"wave") & strcmp(targetLabels,"bottom friction"));
if ~isempty(idx)
    boxDiagram.arrows(idx).intermediatePoints = [16, 5.5+.75; 16, -2; 4.5/2, -2];
    boxDiagram.arrows(idx).LabelPosition = [8 -1.90];
end
% Arrow: M2 tide to damped wave
idx = find(strcmp(sourceLabels,"M2 tide") & strcmp(targetLabels,"damped"+newline+"wave"));
if ~isempty(idx)
    sourceY = sourcesBoxes(end).Position(2)+sourcesBoxes(end).Size(2)/2;
    targetY = dampedBoxes(end).Position(2)+dampedBoxes(end).Size(2)/2;
    sourceX = sourcesBoxes(end).Position(1)+sourcesBoxes(end).Size(1);
    boxDiagram.arrows(idx).intermediatePoints = [sourceX+1, sourceY; sourceX+1, targetY];
    boxDiagram.arrows(idx).LabelPosition = [sourceX+1 sourceY-2];
end    
```