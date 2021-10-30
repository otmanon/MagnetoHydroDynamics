% This script implements the "Fractal Dimension of Dielectric Breakdown",
% the original DBM model introduced by N,P and W in 1984. Unphysical
% because it relies on eta parameter to make nice looking lightning bolts,
% but it is the baseline for all other models basically.
clear; clc;
res = 50;
eta = 1;
pick_max = true;
AI = sub2ind([res, res], floor(0.75*res), floor(0.75*res));        %%indeces of aggregate
SI = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source


[V, F] = create_regular_grid(res);
[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI));
BC(size(BC,1)) = 1;


L = -cotmatrix(V, F);

T = min_quad_with_fixed(L, [], BI, BC);
fig = figure('Position', [0, 0, 1200, 1000]);
clf;
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
s = scatter(V(AI, 1), V(AI, 2), 'filled', 'MarkerFaceColor', 'red');
title("DBM eta 1");
colorbar
colormap(parula(9));
i = 0;
while(true)
    %Get Candidate neighbors
    [AX, AY] = ind2sub(res, AI);
    [N, NI] = getNeighbors([AX, AY], Occ);
    
    
    if (sum(NI == SI) > 0)
        break;
    end
    
    W = T(NI);
    W = W .^ eta;
    W = W ./ sum(W);
    %Sample one of them
    RI = weightedSampling(W, 1);
    if pick_max == true
        [~, RI] = max(W);
    end
    % Set new BC
    AI = [AI; NI(RI)];
    Occ(AI) = 1;
    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;
    %
    T = min_quad_with_fixed(L, [], BI, BC);
    
    t.CData = T;
    s.XData = V(AI, 1);
    s.YData = V(AI, 2);
    
    drawnow;
    if ( mod(i, 3) == 0)
        figgif("./results/DBMDiagMax.gif");
    end
    i = i + 1;
    
end


hold off