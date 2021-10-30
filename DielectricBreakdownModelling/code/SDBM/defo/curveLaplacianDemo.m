clc;
n = 11;
V = zeros(n, 2);
V(:, 1) = 0:1/(n-1):1;
E = zeros(n-1, 2);
E(:, 1) = 1:(n-1);
E(:, 2) = 2:(n);

V2(:, 2) = 0.1:1/(n-1):1.1;
V2(:, 1) = 0.4
E2 = E;

E = [E;7 n+1; E2 + length(V)];
V = [V; V2];


L = curveLaplacian(V, E);

bI = [1; n; length(V)];
bV = [0; 1; 1];

T = min_quad_with_fixed(L, [], bI, bV );

hold on;
p = plot_edges(V, E, 'Color', [0 0 0]);
s = scatter3(V(:, 1), V(:, 2), T, 'filled', 'CData', T);