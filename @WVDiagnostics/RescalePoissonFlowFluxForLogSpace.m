function [logX,logY,Uprime,Vprime] = RescalePoissonFlowFluxForLogSpace(wvd,X,Y,U,V,options)
arguments
    wvd
    X
    Y
    U
    V
    options.shouldOnlyRescaleDirection logical = true
end
logX = log10(X);
logY = log10(Y);
if options.shouldOnlyRescaleDirection == true
    r = 1;
else
    r = sqrt( (1./log10(X)).^2 + (1./log10(Y)).^2 );
end
theta = atan2(abs(log10(X)),abs(log10(Y)));
Uprime = U.*r.*cos(theta);
Vprime = V.*r.*sin(theta);
end