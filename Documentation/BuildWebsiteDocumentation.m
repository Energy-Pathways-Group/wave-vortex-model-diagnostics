cd(fileparts(which(mfilename)))

buildFolder = '../docs/';
sourceFolder = './WebsiteDocumentation/';

copyfile(sourceFolder,buildFolder);

websiteRootURL = "wave-vortex-model-diagnostics/";

classFolderName = 'Class documentation';
websiteFolder = "classes";
classDoc = ClassDocumentation('WVDiagnostics',websiteRootURL=websiteRootURL,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=1);
classDoc.writeToFile();