function [enstrophy, t] = exactEnstrophyOverTime(self, options)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter = @(v) v;
else
    filter = @(v) v(options.timeIndices);
end
enstrophy = filter(self.diagfile.readVariables('enstrophy_apv'));
t = filter(self.diagfile.readVariables('t'));
end