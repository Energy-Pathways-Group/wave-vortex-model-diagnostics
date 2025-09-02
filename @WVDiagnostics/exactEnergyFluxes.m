function energy_fluxes = exactEnergyFluxes(self)
% Return the exact energy flux from the forcing terms, [j kRadial t]
%
% Returns the exact energy fluxes from external forcing and nonlinear
% advection.
%
% - Topic: Fluxes over time, [t 1]
% - Declaration: forcing_fluxes = exactEnergyFluxes(self)
% - Returns forcing_fluxes: struct array with exact fluxes
arguments
    self WVDiagnostics
end
forcingNames = self.wvt.forcingNames;
nForcings = length(forcingNames);
if self.diagfile.hasVariableWithName("E_antialias_filter")
    nForcings = nForcings + 1;
end
energy_fluxes(nForcings) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    energy_fluxes(iForce).name = name;
    energy_fluxes(iForce).fancyName = forcingNames(iForce);
    energy_fluxes(iForce).te = self.diagfile.readVariables("E_" + name);
end

if self.diagfile.hasVariableWithName("E_antialias_filter")
    energy_fluxes(iForce+1).name = "antialias_filter";
    energy_fluxes(iForce+1).fancyName = "antialias filter";
    energy_fluxes(iForce+1).te = self.diagfile.readVariables("E_antialias_filter");
end
end