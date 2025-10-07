function [M_wwg, F_wwg, ks, js] = quadraticEnergyMirrorTriadFluxes2D(self,options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf
    options.mirrorTriad {mustBeMember(options.mirrorTriad, ["wwg","ggw"])} = "wwg";
end

groupName = "mirror-flux-2d-" + options.mirrorTriad;
tName = "t_" + options.mirrorTriad;
if options.mirrorTriad=="wwg"
    triadPrimaryName = "F_wwg_js_ks";
    triadMirrorName = "pi_w_wwg_js_ks";
else
    triadPrimaryName = "F_ggw_js_ks";
    triadMirrorName = "pi_g_ggw_js_ks";
end

if ~self.diagfile.hasGroupWithName(groupName)
    error("Unable to find the group "+groupName+". You may need to call -create2DMirrorFluxes first.");
end
group = self.diagfile.groupWithName(groupName);

js = group.readVariables("js");
ks = group.readVariables("ks");
F = group.readVariables(triadPrimaryName);
Pi = group.readVariables(triadMirrorName);

if ~isinf(options.timeIndices)
    t_diag = self.t_diag;
    t_diag = t_diag(options.timeIndices);
    t = group.readVariables(tName);
    [found,timeIndices] = ismember(t_diag,t);
    if any(~found)
        error("Some of the time indices you requested were not find. You requested " + length(t_diag) + " indices, but only " + sum(found) + " were found.");
    end
else
    timeIndices=1:size(Pi,3);
end

F_wwg = F(:,:,timeIndices);
Pi = Pi(:,:,timeIndices);

Pi_big = zeros([size(Pi,1)+1 size(Pi,2)+1 size(Pi,3)]);
Pi_big(2:end,2:end,:) = Pi;

M_wwg = diff(diff(Pi_big,1,1),1,2);


end