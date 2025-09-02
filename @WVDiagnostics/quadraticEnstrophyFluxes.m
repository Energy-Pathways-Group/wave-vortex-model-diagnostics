function enstrophy_fluxes = quadraticEnstrophyFluxes(self)
% Return the enstrophy flux from the forcing terms
%
% Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
%
% - Topic: Core function — spatial temporal
% - Declaration: forcing_fluxes = quadraticEnergyFluxes(options)
% - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total].
% - Returns forcing_fluxes: an array of structs
arguments
    self WVDiagnostics
end
forcingNames = self.wvt.forcingNames;
nForcings = length(forcingNames);
if self.diagfile.hasVariableWithName("Z0_antialias_filter")
    nForcings = nForcings + 1;
end
enstrophy_fluxes(nForcings) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    enstrophy_fluxes(iForce).name = name;
    enstrophy_fluxes(iForce).fancyName = forcingNames(iForce);
    enstrophy_fluxes(iForce).Z0 = self.diagfile.readVariables("Z0_" + name);
end

if self.diagfile.hasVariableWithName("Z0_antialias_filter")
    enstrophy_fluxes(iForce+1).name = "antialias_filter";
    enstrophy_fluxes(iForce+1).fancyName = "antialias filter";
    enstrophy_fluxes(iForce+1).Z0 = self.diagfile.readVariables("Z0_antialias_filter");
end
end