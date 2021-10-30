function [I, J] = world2grid(VX, VY, size)
%WORLD2GRID Convers list of vertices V to list of indices i, j where i is
%column index and j is row index. N is the number of entries per row/col.
    I = (floor((VX ./ size )-0.5*size)) + 1;
    J = (floor((VY ./ size )-0.5*size)) + 1;
end

