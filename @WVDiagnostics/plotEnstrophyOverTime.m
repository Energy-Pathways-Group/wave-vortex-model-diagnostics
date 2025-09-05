function fig = plotEnstrophyOverTime(self,options)
% Plot enstrophy over time
%
% Plots the quadratic and APV enstrophy as a function of time.
%
% - Topic: Figures (over time)
% - Declaration: fig = plotEnstrophyOverTime(self,options)
% - Parameter options.visible: figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.visible = "on"
    options.timeIndices = Inf;
end
[Z_quadratic, t] = self.quadraticEnstrophyOverTime(timeIndices=options.timeIndices);
[Z_apv, ~] = self.exactEnstrophyOverTime(timeIndices=options.timeIndices);

fig = figure(Visible=options.visible);
plot(t/self.tscale,Z_quadratic/self.zscale,LineWidth=2), hold on
plot(t/self.tscale,Z_apv/self.zscale,LineWidth=2)
legend('quadratic','apv')

xlabel("time (" + self.tscale_units + ")")
ylabel("enstrophy (" + self.zscale_units + ")")
xlim([min(t) max(t)]/self.tscale);
end