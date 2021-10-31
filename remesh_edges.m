function [V, E] = remesh_edges(V, E, max_length, min_length)
%REMESH_EDGES Summary of this function goes here
%   Detailed explanation goes here
    lengths = edge_lengths(V, E);
    
    
    too_long_ind = lengths > max_length;
    while(sum(too_long_ind))
        midP = 0.5*( V(E(too_long_ind, 1), :) + V(E(too_long_ind, 2), :));
        new_indices = size(V, 1) + (1:size(midP, 1));
        newE = E(too_long_ind, :);
        E(too_long_ind, 1) = new_indices;
        newE(:, 2) = new_indices; 
        E = [E; newE];
        V = [V; midP];

        lengths = edge_lengths(V, E);
        too_long_ind = lengths > max_length;
    end
    
   
    
end
