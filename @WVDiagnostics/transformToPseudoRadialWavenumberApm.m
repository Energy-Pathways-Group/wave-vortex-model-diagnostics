function [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)     
% transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial
%
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

wvt = self.wvt;
if size(wvt.h_pm,2) > 1
    kj2 = 1./self.Lr2_pm;
    kj2(1,:) = 0;
    kr2 = repmat(reshape(self.kRadial.^2,1,[]),[length(self.j) 1]);
    Kh = sqrt(kj2 + kr2);
else
    [kj,kr] = ndgrid(self.jWavenumber,self.kRadial);
    Kh = sqrt(kj.^2 + kr.^2);
end

k = self.kPseudoRadial;
dk = k(2)-k(1);
nK = length(k);

varargout = cell(size(varargin));
spectralMatrixSize = size(Kh);
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= spectralMatrixSize(2)
        error('The input matrix must be of size [Nj NkRadial]');
    end
    
    varargout{iVar} = zeros([nK 1]);
end

totalIndices = false(size(Kh));
for iK = 1:1:nK
    indicesForK = k(iK)-dk/2 <= Kh & Kh < k(iK)+dk/2;
    for iVar=1:length(varargin)
        varargout{iVar}(iK) =  sum(varargin{iVar}(indicesForK));
    end
    totalIndices = totalIndices | indicesForK;
end

end