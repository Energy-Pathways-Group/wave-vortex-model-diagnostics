function BuildWebsiteDocumentation(options)
arguments
    options.rootDir = ".."
end
buildFolder = fullfile(options.rootDir,"docs");
sourceFolder = fullfile(options.rootDir,"Documentation","WebsiteDocumentation");

copyfile(sourceFolder,buildFolder);

websiteRootURL = "wave-vortex-model-diagnostics/";

classFolderName = 'Class documentation';
websiteFolder = "classes";
classDoc = ClassDocumentation('WVDiagnostics',websiteRootURL=websiteRootURL,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=1);
classDoc.writeToFile();
end