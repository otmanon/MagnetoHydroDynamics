function [I, J] = world2grid(V, size)
%WORLD2GRID Convers list of vertices V to list of indices i, j where i is
%column index and j is row index. N is the number of entries per row/col.
    I = (floor(V(:, 1) ./ size - 0.5)) + 1;
    J = (floor(V(:, 2) ./ size  - 0.5)) + 1;
end

