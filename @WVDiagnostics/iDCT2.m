function matrix = iDCT2(N)
% InverseCosineTransformMatrix_DCT2  Inverse of forward DCT-II matrix
%
% This inverts CosineTransformForwardMatrix_DCT2 when the forward uses the
% 2/N scaling above (with no extra row halving).
%
% Implementation: DCT-III with DC (first) COLUMN halved.

matrix = zeros(N,N);

for j = 1:N
    for k = 1:N
        matrix(j,k) = cos(pi*(j-0.5)*(k-1)/N);
    end
end

% Halve the DC column
matrix(:,1) = matrix(:,1)/2;

return
