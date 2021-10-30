function [GX, GY, edgeGX, edgeGY] = macGradient(f,h)
%FDGRADIENT Calculates gradient of f using a staggered MAC grid
%   GX, GY -  the gradient at the centers of the cells of f
%   edgeGX, edgeGY - the gradients on the edges of the staggered grid of f
    res = size(f, 1);
    edgeGX = zeros(res, res+1); %vertical edges associated with GX have one more column
    edgeGY = zeros(res + 1, res);
    
    edgeGX(:, [1,end]) = 0; %0-neumann BC. 
    edgeGY([1, end], :) = 0;
    
    edgeGX(:, 2:end-1) = (f(:, 2:end) - f(:, 1:end-1))/h; %assuming positive x dir increases with col
    edgeGY(2:end-1, :) = (f(2:end, :) - f(1:end-1, :))/h;
    
    GX = zeros(size(f));
    GY = zeros(size(f));
    
    GX(:, :) = 0.5* (edgeGX(:, 2:end) + edgeGX(:, 1:end-1));
    GY(:, :) = 0.5* (edgeGY(2:end, :) + edgeGY(1:end-1, :));
end

