clear; clc;
res = 100;

[V, F] = create_regular_grid(res);
AIo = sub2ind([res, res], floor(res), floor(res));        %%indeces of aggregate
SIo = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source

eta = 10;
pick_max = true;
AIDBM = DBM(res, eta, AIo, SIo, pick_max);

mu = 1;
sig = 0.0;
AIDAM = DAM(res, mu, sig, AIo, SIo);
 

fig = figure('Position', [0, 0, 2000, 1000]);
clf;
s1 = subplot(1, 2, 1);
hold on;
scatter(V(AIDBM, 1), V(AIDBM, 2));
title("DBM");
colorbar

s2 = subplot(1, 2, 2);
hold on;
scatter(V(AIDAM, 1), V(AIDAM, 2));
title("DAM");
colorbar

sgtitle('Average Occupancy for eta=10, over 100 Growths');