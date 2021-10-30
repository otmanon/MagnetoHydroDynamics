clc; clf;
clear; clc;

[V, E] = readOBJ("sampleLightning2D.obj");
V = V(:, 1:2);
E = E(:, 1:2);

bI = boundary_faces(E);
bV = V(bI, 1);

M = massmatrix(V, E, 'voronoi');
delem = diag(M);
invDelem = 1./delem;
invM = diag(delem);
invM = kron(invM, eye(2));
L = curveLaplacian(V, E);
L = kron(L, eye(2));

hold on;
p = plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
s = scatter(V([1; length(V)], 1), V([1; length(V)], 2), 'filled');
xlim([10 120]); ylim([0 70]);
colorbar;
caxis([-1 1]);
disp = 1;

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

