cd(fileparts(which(mfilename)))

buildFolder = '../docs/';
sourceFolder = './WebsiteDocumentation/';

copyfile(sourceFolder,buildFolder);

classFolderName = 'Class documentation';
websiteFolder = 'classes';
classDoc = ClassDocumentation('WVDiagnostics',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=1);
classDoc.writeToFile();