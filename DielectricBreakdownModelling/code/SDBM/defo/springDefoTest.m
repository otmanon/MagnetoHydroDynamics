clc; clf;
clear; clc;

n = 11;
V = zeros(n, 2);
V(:, 1) = 0:1/(n-1):1;
E = zeros(n-1, 2);
E(:, 1) = 1:(n-1);
E(:, 2) = 2:(n);
Vorig = V;
lengths = edge_lengths(Vorig, E);

hold on;
xlim([-0.5 1.5]); ylim([-1 1]);
plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
scatter(V(:, 1), V(:, 2), 'filled');


bI = [1; length(V)];
bV = [V(1, :); V(length(V), :) + [-0.5 -0.5] ];

M = massmatrix(V, E, 'voronoi');
delem = diag(M);
invDelem = 1./delem;
invM = diag(delem);
L = curveLaplacian(V, E);

k = 10000000; % stiffness
%Graph Laplacian
A = adjacency_matrix(E)*k;
D = full(sum(A, 2)); %degree matrix
G = diag(D) - A;

%Now build At.
edgeInd = 1:length(E);
ijv = [E(:, 1) edgeInd' ones(length(E), 1);
        E(:, 2) edgeInd' -ones(length(E), 1)];
AT = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), length(V), length(E));
J = k*AT;
%build graph laplacian
%G =  k*J*J';

%energy is given by: 0.5 x'Gx + x'L'invM L x - x'L'invM L x0 - xTJd
% need to do local/global optimization.

f = repelem([0 -1], length(V), 1);

%compute current d
precompute = [];
for i=1:10000
    d = V(E(:, 1), :) - V(E(:, 2), :);
    d = d ./ vecnorm(d, 2, 2);
    d = lengths.*d;

    [x, precompute] = min_quad_with_fixed(0.5* G , -J * d + f , bI, bV, [], [], precompute );
    V = x;
end;

hold on;
xlim([-0.5 1.5]); ylim([-1 1]);
p = plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
s = scatter(V(:, 1), V(:, 2), 'filled');

edge_lengths(V, E) - lengths