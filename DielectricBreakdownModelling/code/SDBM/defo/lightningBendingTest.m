clc; clf;
clear; clc;

[V, E] = readOBJ("sampleLightning2D.obj");
V = V(:, 1:2);
Vorig = V;
E = E(:, 1:2);

lengths = edge_lengths(Vorig, E);

hold on;
xlim([10 80]); ylim([-10 80]);
plot_edges(V, E, 'Color', [1 0 0], 'LineWidth', 2);
scatter(V(:, 1), V(:, 2), 'filled', 'red');


bI = [1; length(V)];
bV = [V(1, :); V(length(V), :) + [10 -10] ];

M = massmatrix(V, E, 'voronoi');
delem = diag(M);
invDelem = 1./delem;
invM = diag(delem);
L = curveLaplacian(V, E);

k = 100000; % stiffness
%Graph Laplacian
A = adjacency_matrix(E)*k;
D = full(sum(A, 2)); %degree matrix
G = diag(D) - A;

%Now build J.
edgeInd = 1:length(E);
ijv = [E(:, 1) edgeInd' k*ones(length(E), 1);
        E(:, 2) edgeInd' -k*ones(length(E), 1)];
J = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), length(V), length(E));

%energy is given by: 0.5 x'Gx + x'L'invM L x - x'L'invM L x0 - xTJd
% need to do local/global optimization.

%compute current d
precompute = [];
for i=1:100
    d = V(E(:, 1), :) - V(E(:, 2), :);
    d = d ./ vecnorm(d, 2, 2);
    d = lengths.*d;

    [x, precompute] = min_quad_with_fixed(0.5*G  + L'*invM*L, -J * d - L'*invM*L*Vorig, bI, bV );
    V = x;
end;

hold on;
axis equal;
xlim([10 80]); ylim([-10 80]);
p = plot_edges(V, E, 'Color', [0 0 1], 'LineWidth', 2);
s = scatter(V(:, 1), V(:, 2), 'filled', 'blue');

final_lengths = edge_lengths(V, E)
