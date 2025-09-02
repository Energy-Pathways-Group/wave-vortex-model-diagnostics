function enstrophy_fluxes = exactEnstrophyFluxes(self)
% Return the available potential enstrophy flux from the forcing terms, [j kRadial t]
%
% Reads from the diagnostics file and returns an array of structs with
% fields name, fancyName, and a field for each energy reservoir with size
% [j kRadial t]. This includes the nonlinear advection term.
%
% - Topic: Core function — spatial temporal
% - Declaration: forcing_fluxes = exactEnstrophyFluxes(options)
% - Returns enstrophy_fluxes: an array of structs
arguments
    self WVDiagnostics
end
forcingNames = self.wvt.forcingNames;
nForcings = length(forcingNames);
if self.diagfile.hasVariableWithName("Z_antialias_filter")
    nForcings = nForcings + 1;
end
enstrophy_fluxes(nForcings) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    enstrophy_fluxes(iForce).name = name;
    enstrophy_fluxes(iForce).fancyName = forcingNames(iForce);
    enstrophy_fluxes(iForce).Z0 = self.diagfile.readVariables("Z_" + name);
end

if self.diagfile.hasVariableWithName("Z_antialias_filter")
    enstrophy_fluxes(iForce+1).name = "antialias_filter";
    enstrophy_fluxes(iForce+1).fancyName = "antialias filter";
    enstrophy_fluxes(iForce+1).Z0 = self.diagfile.readVariables("Z_antialias_filter");
end
end