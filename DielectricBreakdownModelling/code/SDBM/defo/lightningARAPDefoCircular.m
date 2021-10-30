clc; clf;
clear; clc;

res = 200;
circleV = circleVerts(res, 10);
[V, E] = readOBJ("sampleLightning2D.obj");
V = V(:, 1:2);
Vorig = V;
E = E(:, 1:2);

lengths = edge_lengths(Vorig, E);


M = massmatrix(V, E, 'voronoi');
delem = diag(M);
invDelem = 1./delem;
invM = diag(delem);
L = curveLaplacian(V, E);

k = 1; % stiffness
%Graph Laplacian
A = adjacency_matrix(E);
D = full(sum(A, 2)); %degree matrix
G = diag(D) - A;

%Now build J.
edgeInd = 1:length(E);
ijv = [E(:, 1) edgeInd' ones(length(E), 1);
        E(:, 2) edgeInd' -ones(length(E), 1)];
J = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), length(V), length(E));

%energy is given by: 0.5 x'Gx + x'L'invM L x - x'L'invM L x0 - xTJd
% need to do local/global optimization.
bI = [1; length(V)];
bV = [V(1, :); V(length(V), :) + [10 -10] ];



f = 0;
% f = repelem([0 1], length(V), 1);

 hold on;
 xlim([15 70]); ylim([-10, 60]);
 pred = plot_edges(V, E, 'Color', [1 0 0], 'LineWidth', 2);
 sred = scatter(V(:, 1), V(:, 2), 'filled', 'red');
 title("Stretching + Bending resistant lightning");
%figgif("lightningAll.gif");
V2 = V;
V2orig = Vorig;
alpha = 0.5;
beta = 1 - alpha;
precompute = [];
for step=1:res
    V = Vorig;
    V2 = V2orig;
    bI = [1; length(V)];
    bV = [V(1, :); V(length(V), :) + circleV(step, :) ];
    
    
    bending_energy_quad = alpha*L'*invM*L;
    bending_energy_linear = -alpha*2*L'*invM*L*Vorig;
    
    stretching_energy_quad = beta*0.5*k*G;
%    
    %compute current d
    for i=1:10
        d = V(E(:, 1), :) - V(E(:, 2), :);
        d = d ./ vecnorm(d, 2, 2);
        d = lengths.*d;
        stretching_energy_linear = -beta*k*J * d;
        [x, precompute] = min_quad_with_fixed(bending_energy_quad  + stretching_energy_quad,...
         stretching_energy_linear + bending_energy_linear, bI, bV, [], [], precompute);
        V = x;
    end;

    if (step > 1)
         delete(pblue); delete(sblue);
    end;

     pblue = plot_edges(V, E, 'Color', [0 0 1], 'LineWidth', 2);
     sblue = scatter(V(:, 1), V(:, 2), 'filled', 'blue');
%   

    drawnow;
    figgif("lightningARAP.gif");
    
end;
    


function V = circleVerts(res, rad);
    theta = 0:2*pi/res:2*pi;
    
    x = cos(theta)';
    y = sin(theta)';
    V = rad.*[x y];
end
