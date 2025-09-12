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
forcingNames = self.forcingNames;
enstrophy_fluxes(length(forcingNames)) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    enstrophy_fluxes(iForce).name = name;
    enstrophy_fluxes(iForce).fancyName = forcingNames(iForce);
    enstrophy_fluxes(iForce).Z0 = self.diagfile.readVariables("Z0_" + name);
end

end