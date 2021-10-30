clear; clc;
res = 25;
eta = 10;
pick_max = false;

[V, F] = create_regular_grid(res);
AIo = sub2ind([res, res], floor(res), floor(res));        %%indeces of aggregate
SIo = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source

Occ1 = zeros(size(V, 1), 1);
Occ2 = zeros(size(V, 1), 1);
for i=1:100
     AI2 = DBM(res, eta, AIo, SIo, pick_max);
     AI = AI2;
     Occ1(AI) = Occ1(AI) + 1;
 end

sampling_density = 10000
for i=1:100
    AI2 = DBMDistanceIsosurface(res, eta, AIo, SIo, pick_max, sampling_density);
    AI = AI2;
    Occ2(AI) = Occ2(AI) + 1;
end
fig = figure('Position', [0, 0, 2000, 1000]);
clf;
s1 = subplot(1, 2, 1);
hold on;
t = tsurf(F, V, 'CData', Occ1 ,  fphong, falpha(1, 0));
title("Sampling Grid Neighbors");
colorbar

s1 = subplot(1, 2, 2);
hold on;
t = tsurf(F, V, 'CData', Occ2 ,  fphong, falpha(1, 0));
title("Sampling Grid Distance Isosurface, s.d. 10000");
colorbar

sgtitle('Average Occupancy for eta=10, over 100 Growths');