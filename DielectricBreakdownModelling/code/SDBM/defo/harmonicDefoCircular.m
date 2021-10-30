clc;clf;
clear; clc;

[V, E] = readOBJ("sampleLightning2D.obj");
V = V(:, 1:2);
Vorig = V;
E = E(:, 1:2);

bI = boundary_faces(E);
bV = V(bI, 1);


L = curveLaplacian(V, E);
L = kron(L, eye(2));

hold on;
xlim([10 120]); ylim([-30 70]);
colorbar;
caxis([-1 1]);
title("harmonic deformation");
figgif("harmonic_defo_circular.gif");
res = 200;
circleV = circleVerts(res, 1);
for step=1:res
    
    bIflat = [1; 2; 2*length(V) - 1; 2*length(V)];
    
    disp = circleV(step, :);
 
    bV = [0; 0; -20*disp'];

    d = min_quad_with_fixed(L, [], bIflat, bV );

    dMat = vec2mat(d, length(V), 2); 
    V = Vorig + dMat;
    if (step > 1)
        delete(p); delete(s);
    end;
    p = plot_edges(V, E, 'Color', [0 0 0], 'LineWidth', 2);
    s = scatter(V([1; length(V)], 1), V([1; length(V)], 2), 'filled');
    drawnow; 
    
    figgif("harmonic_defo_circular.gif");
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


function V = circleVerts(res, rad);
    theta = 0:2*pi/res:2*pi;
    
    x = cos(theta)';
    y = sin(theta)';
    V = [x y];
end
