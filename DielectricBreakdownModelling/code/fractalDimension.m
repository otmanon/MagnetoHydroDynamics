clear; clc;
res = 200;
eta = 1;
pick_max = false;

[V, F] = create_regular_grid(res);
AIo = sub2ind([res, res], floor(res), floor(res));        %%indeces of aggregate
SIo = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source

Occ1 = zeros(size(V, 1), 1);


% sampling_density = 10000
% for i=1:100
%     AI2 = DBMDistanceIsosurface(res, eta, AIo, SIo, pick_max, sampling_density);
%     AI = AI2;
%     Occ2(AI) = Occ2(AI) + 1;
% end

AI = DBM(res, eta, AIo, SIo, pick_max);
 

 radii = [1:4:res-10] .* 1/res;
 Nr = zeros(1, length(radii));
for i=1:length(radii)
    r = radii(i);
    idx = rangesearch(V(AI, :), V(AI, :), r);
    numNeigh = cellfun(@length, idx);
    avgNumNeigh = mean(numNeigh);
    Nr(i) = avgNumNeigh;
    
end
lognr = log(Nr);
loginvr = log(1./radii);
p= polyfit( loginvr, lognr,1);
x = min(loginvr):(max(loginvr) -min(loginvr)) /10 :max(loginvr);
y = p(1) * x + p(2);

fig = figure('Position', [0, 0, 1200, 1000]);
clf;
s1 = subplot(1, 2, 1);
hold on;
t = tsurf(F, V, 'CData', Occ1 ,  fphong, falpha(1, 0));
s = scatter(V(AI, 1), V(AI, 2), 'filled',  'MarkerFaceColor', 'red');
title("Aggregate");
colorbar

s2 = subplot(1, 2, 2);

hold on;
plot( loginvr, lognr,  'b');
plot(x, y, 'r');
xtitle("log(1/r)");
ytitle("log(N(r))");
title("Fractal Dimension");

