function [VX, VY, nodeVX, nodeVY] = regular_grid(res, size)
    x = 0:size:size*(res-1);
    y = x;
    [VX, VY] = meshgrid(x, y);
    VX = VX + 0.5*size;
    VY = VY + 0.5*size;
    %Again for the nodes
    x = 0:size:size*(res);
    y = x;
    [nodeVX, nodeVY] = meshgrid(x, y);

end