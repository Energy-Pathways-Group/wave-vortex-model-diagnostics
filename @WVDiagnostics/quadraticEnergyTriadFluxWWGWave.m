function wwg = quadraticEnergyTriadFluxWWGWave(self)
[fpath,fname,~] = fileparts(self.wvpath);
if isempty(fpath)
    fpath = pwd;
end
wwgTriadPath = fullfile(fpath,strcat(fname,"-wwg-aa.nc"));
if ~exist(wwgTriadPath,"file")
    wwgTriadPath = fullfile(fpath,strcat(fname,"-wwg.nc"));
    if ~exist(wwgTriadPath,"file")
        error("Unable to find the wwg file. You can create one using the -createWWGTriadDiagnostic function.");
    end
end
ncfile = NetCDFFile(wwgTriadPath,shouldReadOnly=true);
t = ncfile.readVariables("t");

wavewave_jkt = ncfile.readVariables("wavewave_jk");
wavewave_flux_t = zeros(size(wavewave_jkt));
wavewave_flux_t(1,1,:) = wavewave_jkt(1,1,:);
wavewave_flux_t(2:end,1,:) = diff(wavewave_jkt(:,1,:),1,1);
wavewave_flux_t(1,2:end,:) = diff(wavewave_jkt(1,:,:),1,2);
wavewave_flux_t(2:end,2:end,:) = diff(diff(wavewave_jkt,1,1),1,2);

% If the wvd had explicit de-aliasing, but this diagnostic did not, then
% this will zeropad.
wwg = zeros(length(self.j),length(self.kRadial),length(t));
wwg(1:size(wavewave_flux_t,1),1:size(wavewave_flux_t,2),:) = wavewave_flux_t;

end