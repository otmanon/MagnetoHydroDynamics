clc;
clear; clc;

[V, E] = readOBJ("sampleLightning2D.obj");
V = V(:, 1:2);
E = E(:, 1:2);

bI = boundary_faces(E);
bV = V(bI, 1);

bI = [1; length(V)];
bV = [0; 1];

L = curveLaplacian(V, E);


T = min_quad_with_fixed(L, [], bI, bV );

hold on;
p = plot_edges(V, E, 'Color', [0 0 0]);
s = scatter3(V(:, 1), V(:, 2), T, 'filled', 'CData', T);