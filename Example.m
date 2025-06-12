rootdir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
figureFolder = rootdir + "figures";
if ~exist(figureFolder, 'dir')
       mkdir(figureFolder)
end

experiment = configureDictionary("string","string");
% experiment("hydrostatic, mean flow") = "run1_icR_iner0_tide0_lat32_geo0065_N0052_hydrostatic_res256.nc";
% experiment("hydrostatic, mean flow + nio + M2") = "run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res256.nc";
% experiment("nonhydrostatic, mean flow + nio + M2") = "run18_icR_iner07_tide014_lat32_geo0065_N0052_boussinesq_res256.nc";
experiment("nonhydrostatic, mean flow + nio + M2 512") = "run99_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512.nc";

%%
custom_names = configureDictionary("string","string");
custom_names("quadratic_bottom_friction") = "bottom friction";
custom_names("vertical_diffusivity") = "diffusivity";
custom_names("adaptive_damping") = "damping";
custom_names("inertial_forcing") = "NIO";
custom_names("M2_tidal_forcing") = "M2 tide";
custom_names("geostrophic_mean_flow") = "mean flow";
custom_names("te_gmda") = "geostrophic";
timeIndices=750:1001;

%%
keys = experiment.keys;
for iExperiment = 1:length(keys)
    experimentName = keys(iExperiment);
    tmp = split(experiment(experimentName),"_");
    experimentNumber = tmp(1);

    experimentFolder = figureFolder + "/" + experimentNumber;
    if ~exist(experimentFolder, 'dir')
        mkdir(experimentFolder)
    end

    wvd = WVDiagnostics(rootdir + experiment(experimentName));

    fig = wvd.plotEnstrophyOverTime(visible="off");
    exportgraphics(fig,experimentFolder + "/" + "enstrophy_vs_time.png")
    
    fig = wvd.plotEnergyForReservoirOverTime(visible="off");
    exportgraphics(fig,experimentFolder + "/" + "energy_vs_time.png")

    % fig = wvd.plotMooringRotarySpectrum(visible="off");
    % exportgraphics(fig,experimentFolder + "/" + "mooring_spectrum.png")

    title = experimentName + " days " + round(min(wvd.t_diag(300:441))/86400) + "-" + round(max(wvd.t_diag(300:441))/86400);
    fig = wvd.plotSourcesSinksReservoirsDiagram(customNames=custom_names,shouldShowUnits=true,timeIndices=300:441,title=title,visible="off");
    exportgraphics(fig,experimentFolder + "/" + "energy_flux_diagram.png",Resolution=300)
end