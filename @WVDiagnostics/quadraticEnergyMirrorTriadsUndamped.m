function [E0_ggw,Epm_wwg] = quadraticEnergyMirrorTriadsUndamped(self,options)
arguments
    self WVDiagnostics
    options.timeIndices
end
[fpath,fname,~] = fileparts(self.wvpath);
if isempty(fpath)
    fpath = pwd;
end
wwgTriadPath = fullfile(fpath,strcat(fname,"-damped-aa.nc"));
if ~exist(wwgTriadPath,"file")
    wwgTriadPath = fullfile(fpath,strcat(fname,"-damped.nc"));
    if ~exist(wwgTriadPath,"file")
        error("Unable to find the damped file. You can create one using the -create.. function.");
    end
end

if exist(wwgTriadPath,"file")
    %%%%%%%%%%%%%%%%%%
    % If existing file
    %%%%%%%%%%%%%%%%%%
    ncfile = NetCDFFile(wwgTriadPath);
    [E0_ggw_t,Epm_wwg_t] = ncfile.readVariables("E0_ggw_nodamp","Epm_wwg_nodamp");
    if ~isinf(options.timeIndices)
        E0_ggw_t = E0_ggw_t(:,:,options.timeIndices);
        Epm_wwg_t = Epm_wwg_t(:,:,options.timeIndices);
    end
    E0_ggw = zeros(length(self.j),length(self.kRadial),size(E0_ggw_t,3));
    E0_ggw(1:size(E0_ggw_t,1),1:size(E0_ggw_t,2),:) = E0_ggw_t;
    Epm_wwg = zeros(length(self.j),length(self.kRadial),size(Epm_wwg_t,3));
    Epm_wwg(1:size(Epm_wwg_t,1),1:size(Epm_wwg_t,2),:) = Epm_wwg_t;
else
    error("unable to find the appropriate file.")
end
end