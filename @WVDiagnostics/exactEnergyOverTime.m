function [energy, t] = exactEnergyOverTime(self, options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter = @(v) v;
else
    filter = @(v) v(options.timeIndices);
end
[ke,ape] =self.diagfile.readVariables('ke','ape');
energy = filter(ke+ape);
t = filter(self.diagfile.readVariables('t'));
end