clc; clear;


res = 100;                  %number of grid vertices
cell_size = 1;                    %size of each cell... always make it so that the domain is the [0,1]x[0, 1] unit square


eps = 0.1;                 %minimum radius for WoS solve.
numIter=1000;

[nodeV, nodeF] = create_regular_grid(res);

%FEM params
C = cotmatrix(nodeV, nodeF);
bE = boundary_faces(nodeF);
bI = unique(bE);

leftI = nodeV(bI, 1) == min(nodeV(:, 1));
bV = ones(length(bI), 1); % boundary value is just x coordinate
bV(leftI) = 0;
TFEM = min_quad_with_fixed(C, [], bI, bV);

%WoS
[V, IM, IN] = remove_unreferenced(nodeV, bE);
E = IM(bE);
bV = ones(length(V), 1);
leftI = V(:, 1) == min(V(:, 1));
bV(leftI) = 0;


[TWoS, gradTWoS] = WoS(nodeV, V, E, bV, 1e-5, numIter );


fig = figure('Position', [100, 100, 1000, 400]);
subplot(2, 2, 1);
hold on;
colormap(parula(100));
colorbar();
title('FEM');
tsurf(nodeF, nodeV, 'CData', TFEM, fphong, falpha(1, 0));

subplot(2, 2, 2);
hold on;
colormap(parula(100));
colorbar();
title('WoS');
tsurf(nodeF, nodeV, 'CData', TWoS, fphong, falpha(1, 0));

index=55;
scatter(nodeV(index, 1), nodeV(index, 2), 'filled');



subplot(2, 2,  4);

quiver(nodeV(:, 1), nodeV(:, 2), gradTWoS(:, 1), gradTWoS(:, 2));
