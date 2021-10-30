clear; clc;
res = 50;

[V, F] = create_regular_grid(res);
AIo = sub2ind([res, res], floor(res), floor(res));        %%indeces of aggregate
SIo = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source

eta =3;
pick_max = false;
sampling_density = 5000;
r = 0.06;
[AIDBMGrid, AVGrid] = DBM(res, eta, AIo, SIo, pick_max );
[AIDBMIso, AVIso] = DBMDistanceIsosurface(res, eta, AIo, SIo, pick_max, sampling_density );
[AVGrid, AVI, AVJ] = remove_duplicate_vertices(AVGrid, 1e-7);



fig = figure('Position', [0, -500, 1500, 1200]);
fontsize = 18;
clf;
s1 = subplot(2, 2, 1);
hold on;
s = scatter(V(AIDBMGrid, 1), V(AIDBMGrid, 2),  'red', 'filled');
title("Grid Sampled Candidates", 'FontSize', fontsize);
xlim([-0.2, 1.2]);
ylim([-0.2, 1.2]);
sgtitle(strcat('Local Neighborhoods Grid Alignment, eta=3, r=', num2str(r)), 'FontSize', fontsize*2);
scale = 0.01;
[indices,dists] = rangesearch(AVGrid,AVGrid, r);
principalDir = zeros(2*size(AVGrid, 1), 2);
for i=1:size(indices)
    P = AVGrid(indices{i}, :);
    meanP = mean(P); % AV(i, :);
    P = P-meanP;
    Cov = P'*P;
    [EV, D] = eig(Cov);
    [d, ind] = sort(diag(D), 'descend');
    D = D(ind, ind);
    EV = EV(:, ind);
    EV = EV ./ vecnorm(EV, 2, 2);
    principalDir(2*i - 1, :) = EV(:, 1);
    principalDir(2*i, :) = -EV(:, 1);
    quiver(AVGrid(i, 1), AVGrid(i, 2), scale*EV(1, 1), scale*EV(2, 1), 'LineWidth', 2, 'color' ,'blue');
end


xaxis = repelem([1, 0], 2*size(AVGrid, 1), 1);
yaxis = repelem([0, 1], 2*size(AVGrid, 1), 1);

costheta = dot(principalDir, xaxis, 2);
sintheta = dot(principalDir, yaxis, 2);
angle = atan2( sintheta, costheta);

nbins = 36;
bins = -2*pi/(nbins):2*pi/nbins:2*pi-2*pi/(nbins);
%bins = wrapTo2Pi(bins);
s2 = subplot(2,2, 3);
rplot  = rose(angle, bins);
hold on;
title("Alignment", 'FontSize',fontsize);




s3 = subplot(2, 2, 2);
hold on;
sc3 = scatter(V(AIDBMIso, 1), V(AIDBMIso, 2),  'red', 'filled');
title("Iso-surface Sampled Candidates", 'FontSize', fontsize);
xlim([-0.2, 1.2]);
ylim([-0.2, 1.2]);
[indices2,dists2] = rangesearch(AVIso,AVIso, r);
principalDir2 = zeros(2*size(AVIso, 1), 2);
for i=1:size(indices2)
    P = AVIso(indices2{i}, :);
    meanP = mean(P); % AV(i, :);
    P = P-meanP;
    Cov = P'*P;
    [EV, D] = eig(Cov);
    [d, ind] = sort(diag(D), 'descend');
    D = D(ind, ind);
    EV = EV(:, ind);
    EV = EV ./ vecnorm(EV, 2, 2);
    principalDir2(2*i - 1, :) = EV(:, 1);
     principalDir2(2*i, :) = -EV(:, 1);
    quiver(AVIso(i, 1), AVIso(i, 2), scale*EV(1, 1), scale*EV(2, 1), 'LineWidth', 2, 'color' ,'blue');
end


xaxis2 = repelem([1, 0], 2*size(AVIso, 1), 1);
yaxis2 = repelem([0, 1], 2*size(AVIso, 1), 1);

costheta2 = dot(principalDir2, xaxis2, 2);
sintheta2 = dot(principalDir2, yaxis2, 2);
angle2 = atan2( sintheta2, costheta2);
s4 = subplot(2,2, 4);
rplot2  = rose(angle2, bins);
hold on;
title("Alignment", 'FontSize',fontsize);
 