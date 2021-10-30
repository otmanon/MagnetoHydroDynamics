% This method implements the "Fast Simulation of Laplacian Growth" method
% by Kim et. al (2007). This method isn't truly solving the same equations
% as the others, the summing of the Green's Functions neglects imposing of
% proper boundary conditions and results in a difficult to control, "Center of Mass" repelled
% growth

clear; clc;
res = 100;
eta = 1;
[V, F] = create_regular_grid(res);
pick_max = true;
AV = [0.75, 0.75];
AE = [1, 1];
sub = snap_points(AV, V);
AI = sub(1);       %%indeces of aggregate
SI = sub2ind([res, res], floor(res/2), floor(res/2));                 %% First vertex is the source

[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI));
BC(size(BC,1)) = 1;


L = -cotmatrix(V, F);

T = sum_kernels(V, V(AI, :), V(SI, :));
fig = figure('Position', [0, 0, 1200, 1000]);
clf;
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
sa = scatter(V(AI, 1), V(AI, 2), 50, 'filled', 'MarkerFaceColor', 'red');
ss = scatter(V(SI, 1), V(SI, 2), 'filled', 'MarkerEdgeColor', 'black');
title("Green's Function Summing, res 100")
%pa = plot_edges(AV, AE,'red', 'LineWidth', 2);
colorbar;
i = 0;
while(true)
    %Get Candidate neighbors
    [AX, AY] = ind2sub(res, AI);
    [N, NI] = getNeighbors([AX, AY], Occ);
    
    if (size(SI, 1) > 0)
        if (sum(NI == SI) > 0)
            break;
        end
    end

    
    W = T(NI);
    W = W .^ eta;
    W = W ./ sum(W);
    %Sample one of them
    RI = weightedSampling(W, 1);
    if pick_max == true
        [~, RI] = max(W);
    end
    %[AV2, AE2] = growAggregate(AV, AE, V(NI(RI), :));
    %AV = AV2; AE = AE2;
    % Set new BC
    AI = [AI; NI(RI)];
    Occ(AI) = 1;
    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;
    %
    iT = sum_kernels(V, V(NI(RI, :), :), []);
    T = T + iT;
    t.CData = T;
    sa.XData = V(AI, 1);
    sa.YData = V(AI, 2);
   % pa = plot_edges(AV, AE, 'red', 'LineWidth', 2);
    drawnow;
    
    if ( mod(i, 3) == 0)
        figgif("./results/GFKernelMaxTres100.gif");
    end
    i = i + 1;
    
end


hold off