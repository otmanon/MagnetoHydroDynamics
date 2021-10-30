function [L] = fdLaplacian(f, h)
%FDLAPLACIAN Summary of this function goes here
%   Calculates the finite difference laplacian of a scalar field f with zero neumann BC. Only works for 2D 
    
    augF = zeros(size(f) + 2); % augmented f matrix. will put in "ghost values" at edges
    
    augF(2:end-1, 2:end-1) = f;
    
    %     %For points on top edge (but not corners), L(0, 2:end-1) = L(2,
     %2:end-1) because of 0-Neumann BC. Assign ghost value as such
     % Therefore all edges of the BC 
    augF(1, 2:end-1) = f(2, :);  
    
    augF(end, 2:end-1) = f(end-1, :);
    
    augF(2:end-1, 1) = f(:, 2);
    
    augF(2:end-1, end) = f(:, end-1);
    
    L = zeros(size(f));
    
    %Use standard 5 point Laplacian stencil on interior of augF:
    L(:, :) = (augF(1:end-2, 2:end-1) + augF(3:end, 2:end-1)...
            + augF(2:end-1, 1:end-2) + augF(2:end-1, 3:end)...
            - 4*augF(2:end-1, 2:end-1)) / (h*h);
    

end

