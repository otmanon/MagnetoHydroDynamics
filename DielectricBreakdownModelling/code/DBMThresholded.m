% This method implements the Zeller and Weismann (1986) modification, which adds a threshold value
% to the breakdown. If the electric field at a point is below that threshold, it's weight is set to zero 
% and the probability of sampling it is 0. While being more physically
% based, this method introduced a parameter that may stop the growth if the
% electric field is too low. It also suffers from grid artefacts that do
% not go away with resolution.

clear; clc;
res = 50;
eta = 10;
[V, F] = create_regular_grid(res);
threshold  = 0.04;
pick_max = true;
AI = sub2ind([res, res], floor(res), floor(res));        %%indeces of aggregate
SI = sub2ind([res, res], 1, floor(res/2));                 %% First vertex is the source


[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI));
BC(size(BC,1)) = 1;


L = -cotmatrix(V, F);

T = min_quad_with_fixed(L, [], BI, BC);
clf;
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
s = scatter(V(AI, 1), V(AI, 2), 'filled', 'MarkerFaceColor', 'red');
colorbar
while(true)
    %Get Candidate neighbors
    [AX, AY] = ind2sub(res, AI);
    [N, NI] = getNeighbors([AX, AY], Occ);
    
    
    if (sum(NI == SI) > 0)
        break;
    end
    
    W = T(NI);
    W(W<threshold) = 0;
    
    
    W = W .^ eta;
    W = W ./ sum(W);
    %Sample one of them
    if (sum(W) > 0)
        RI = weightedSampling(W, 1);
    else
        continue
    end
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
    
    
end


hold off