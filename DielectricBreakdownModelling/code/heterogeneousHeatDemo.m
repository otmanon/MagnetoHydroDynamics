% This script implements the "Phase-Field MOdel of Dielectric Breakdown in Solids",.
clear; clc;
res = 50;

[V, F] = create_regular_grid(res);



T = zeros(size(V, 1), 1);
stind = sub2ind([res res], 5, 10);
T(V(:, 1)> 0.45 & V(:, 1)< 0.55 & V(:, 2)> 0.45 & V(:, 2)< 0.55) = 1;

s1 = subplot(1, 1, 1);
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
title("Heterogeneous Heat Diffusion Solve");
colorbar();
caxis([0, 1]);
colormap(parula(30));
deltaT = 0.1;

alphas = V(:, 1).^2 .* V(:, 2).^2;

   
L = fdLaplacianMatrix(res, alphas);% cotmatrix(V, F);
step = 0;
while(true)

   laplacian = L*T;
    T = T + deltaT * laplacian;
    T(V(:, 1)> 0.45 & V(:, 1)< 0.55 & V(:, 2)> 0.45 & V(:, 2)< 0.55) = 1;
    t.CData = T;
   
    drawnow;
    if (mod(step, 100) == 0)
        figgif("results/hetHeatSolve3.gif");
    end
    step = step + 1;
end

s2 = subplot(1, 2, 2);
hold on;
colorbar();
%L = cotmatrix(V, F);
bI = unique(boundary_faces(F));
bT = V(bI, 1);   
T = min_quad_with_fixed(L, [], bI, bT);
               % boundary conditions according ot x coord


t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));



hold off