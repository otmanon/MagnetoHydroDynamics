% This script implements the "Phase-Field MOdel of Dielectric Breakdown in Solids",.
clear; clc;
res = 50;
h = 1/(res-1);


[V, F] = create_regular_grid(res);
[Ae, NBC] = neumannGridBC(res, h, 0.25); 


xi = (ceil(0.25*res):floor(0.75*res))';
%yi = (floor(1):floor(1))';
yi = repelem(ceil(res/2), length(xi), 1);
ii = sub2ind([res res], yi, xi);        %starting index of breakdown
phi = zeros(size(V, 1), 1); %potential
s = ones(size(V, 1), 1);    %phase field variable
s(ii) = 0;


botI = find(V(:, 2) < 0.01);
topI = find(V(:, 2) > 0.99);
% BI = [botI; topI; ii];
% BPhi = [20*ones(length(botI), 1); zeros(length(topI), 1); zeros(length(ii), 1)] ;   
%  BI = [res*res/2 - floor(res/4);res*res/2 - floor(0.75*res);  ii];
%  BPhi = [1; -1;  zeros(length(ii), 1)] ;   
 BI = [res*res/2 - floor(res/4);  ii];
 BPhi = [1;   zeros(length(ii), 1)] ;   

eps = ones(length(V), 1);

LPhi = fdLaplacianMatrix(res, eps, h);% cotmatrix(V, F);
phi = min_quad_with_fixed(LPhi, [], BI, BPhi, Ae, NBC);  %solve for potential field;
     
[EX, EY] = macGradient(reshape(phi, [res, res]), 1);

E = [EX(:), EY(:)];
energy = dot(E, E, 2);

fontsize = 24;
fig = figure('Position', [0, 0, 1500, 600]);
s1 = subplot(2, 3, 1);
hold on;
tPlot = tsurf(F, V, 'CData', phi, falpha(1, 0), fphong);
title("Potential DBC", 'FontSize', fontsize);
colorbar();
colormap(parula(100));

s2 = subplot(2, 3, 4);
hold on;
colorbar();
esPlot = tsurf(F, V, 'CData', energy, fphong, falpha(1, 0));
title(strcat("Dirichlet Energy : ", num2str(sum(energy))), 'FontSize', fontsize);


eps = ones(length(V), 1);
eps(ii) = 10;
% BI = [botI; topI];
% BPhi = [20*ones(length(botI), 1); zeros(length(topI), 1)] ;   
BI = [res*res/2 - floor(res/4)];
BPhi = [1] ;   
LPhi = fdLaplacianMatrix(res, eps, h);
phi = min_quad_with_fixed(LPhi, [], BI, BPhi, Ae, NBC);  %solve for potential field;
     
[EX, EY] = macGradient(reshape(phi, [res, res]), 1);
E = [EX(:), EY(:)];
energy = dot(E, E, 2);



s3 = subplot(2, 3, 2);
hold on;
tPlot = tsurf(F, V, 'CData', phi, falpha(1, 0), fphong);
title("Potential Eps", 'FontSize', fontsize);
colorbar();

s4 = subplot(2, 3, 5);
hold on;
colorbar();
esPlot = tsurf(F, V, 'CData', energy, fphong, falpha(1, 0));
title(strcat("Dirichlet Energy : ", num2str(sum(energy))), 'FontSize', fontsize);


eps = ones(length(V), 1);
% BI = [botI; topI];
% BPhi = [20*ones(length(botI), 1); zeros(length(topI), 1)] ;   

 BI = [res*res/2 - floor(res/4)];
 BPhi = [1] ;   
LPhi = fdLaplacianMatrix(res, eps, h);
phi = min_quad_with_fixed(LPhi, [], BI, BPhi, Ae, NBC);  %solve for potential field;
     
[EX, EY] = macGradient(reshape(phi, [res, res]), 1);
E = [EX(:), EY(:)];
energy = dot(E, E, 2);



s5 = subplot(2, 3, 3);
hold on;
tPlot = tsurf(F, V, 'CData', phi, falpha(1, 0), fphong);
title("Potential None", 'FontSize', fontsize);
colorbar();

s6 = subplot(2, 3, 6);
hold on;
colorbar();
esPlot = tsurf(F, V, 'CData', energy, fphong, falpha(1, 0));
title(strcat("Dirichlet Energy : ", num2str(sum(energy))), 'FontSize', fontsize);