clc; clf;
clear; clc;

n = 11;
V = zeros(n, 2);
V(:, 1) = 0:1/(n-1):1;
E = zeros(n-1, 2);
E(:, 1) = 1:(n-1);
E(:, 2) = 2:(n);

Vorig = V;


hold on;
xlim([-0.5 1.5]); ylim([-1 1]);
plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
scatter(V(:, 1), V(:, 2), 'filled');



bI = [1; 2; length(V)];
bV = [V(1, 1) V(1, 2); V(2, 1) V(2, 2); V(length(V), 1)-0.55 V(length(V), 2)-0.55];


M = massmatrix(V, E, 'voronoi');
delem = diag(M);
invDelem = 1./delem;
invM = diag(delem);
L = curveLaplacian(V, E);

x0 = Vorig';
x0 = x0(:);

x = min_quad_with_fixed(L'*invM*L, L'*invM*L*Vorig, bI, bV );
V = x

hold on;
xlim([-0.5 1.5]); ylim([-1 1]);
p = plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
s = scatter(V(:, 1), V(:, 2), 'filled');

colorbar;
caxis([-1 1]);
disp = 1;
drawnow;
figgif("biharmonic_defo_linear.gif");
for step=1:100
    
    bIflat = [1; 2; 2*length(V) - 1; 2*length(V)];
    
  
    if (mod(step, 50) == 0);
        disp = -1*disp;
    end;
    bV = [0; 0; disp; 0];

    d = min_quad_with_fixed(L, [], bIflat, bV );

    dMat = vec2mat(d, length(V), 2); 
    V = V + dMat;
    delete(p); delete(s);
    p = plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
    s = scatter(V([1; length(V)], 1), V([1; length(V)], 2), 'filled', 'CData', dMat([1; length(V)]));
    drawnow; 
    
    figgif("biharmonic_defo_linear.gif");
end

%s = scatter3(V(:, 1), V(:, 2), T, 'filled', 'CData', T);

function mat = vec2mat(vec, rows, cols)
    %assumes vec is organized in row major order
    ind = 1:length(vec);
    ind = mod(ind-1, cols) + 1;
    newVec = [];
    for i=1:cols
        newVec = [newVec; vec(ind==i)];
    end;
    mat = reshape(newVec, rows, cols);
end

