function [div] = macDivergence(edgeGX, edgeGY, h)
%MACDIVERGENCE Calculates divergence of vector field given by staggered MAC
% grid.
%   Detailed explanation goes here
    res = min(size(edgeGX));
    div = zeros(res, res);
    
    %divergence is simply GX2-GX1 + GY2-GY1 / h
    div(:, :) = edgeGX(:, 2:end) - edgeGX(:, 1:end-1) ...
            + edgeGY(2:end, :) - edgeGY(1:end-1, :);
    div = div/h;
        
    
    
end

