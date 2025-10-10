function [M_ggw, kePeAxis] = quadraticEnergyMirrorTriadFluxes1D_kepe(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

val= self.diagfile.readVariables('pi_g_ggw_kepe');
M_ggw= cat(1,zeros(1,size(val,2)),diff(val));

kePeAxis =  reshape(self.diagfile.readVariables('kePeAxis'),[],1);

if isinf(options.timeIndices)
    M_ggw = mean(M_ggw,2);
else
    M_ggw = mean(M_ggw(:,options.timeIndices),2);
end

end