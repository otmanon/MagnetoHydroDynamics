function L = curveLaplacian(V, E)
%CURVELAPLACIAN Summary of this function goes here
%   Detailed explanation goes here
    
    l = edge_lengths(V, E);
    I1 = E(:,1 );
    I2 = E(:, 2);
    
    invl = 1./l;
    ijv = [I1 I1 -invl; I2 I2 -invl; I1 I2 invl; I2 I1 invl ];
    
    
    L = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), length(V), length(V));
    
    
end

