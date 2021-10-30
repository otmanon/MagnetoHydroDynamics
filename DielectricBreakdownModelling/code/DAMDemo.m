% This script implements the "Fractal Dimension of Dielectric Breakdown",
% the original DBM model introduced by N,P and W in 1984. Unphysical
% because it relies on eta parameter to make nice looking lightning bolts,
% but it is the baseline for all other models basically.
clear; clc;
res = 50;
eta = 1;
pick_max = true;
AI = sub2ind([res, res], floor(res*.75), floor(res*.75));        %%indeces of aggregate
SI = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source


[V, F] = create_regular_grid(res);
[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
breakdownField = zeros(size(V, 1), 1);
breakdownField(AI) = 1;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI));
BC(size(BC,1)) = 1;


L = -cotmatrix(V, F);

T = min_quad_with_fixed(L, [], BI, BC);
fig = figure('Position', [0, 0, 1200, 1000]);
clf;
s1 = subplot(1, 1, 1);
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
s1 = scatter(V(AI, 1), V(AI, 2), 'filled', 'MarkerFaceColor', 'red');
title("Potential Field");
colorbar();
colormap(parula(9));

% s2 = subplot(1, 2, 2);
% hold on;
% b = tsurf(F, V, 'CData', breakdownField, fphong, falpha(1, 0));
% s2 = scatter(V(AI, 1), V(AI, 2), 'filled', 'MarkerFaceColor', 'red');
% title("Breakdown Field");
% colorbar();
% colormap(parula(9));

i = 0;
h = 1/res;
lambda = 100;
Lb = 1;
c = 1;
prevAIsize = 0;


noise_mu = 1;
noise_sig = 1.0;
perturbation = (normrnd(noise_mu, noise_sig, size(V, 1), 1)) ;
T = perturbation.*T;
while(true)
    %Get Candidate neighbors

    [AX, AY] = ind2sub(res, AI);
    [N, NI] = getNeighbors([AX, AY], Occ);
    
    if (sum(NI == SI) > 0)
        break;
    end

    alpha = 1/lambda * exp(- 1./(lambda*T(NI)));
    fq = c*(exp(alpha.* Lb) - 1);
    %for all neighbors find the first one that gets to 1
    diff = ones(size(NI, 1), 1) - breakdownField(NI);
   
    stepsNeeded = diff ./ fq;
    stepsNeeded(diff == 0) = 0;
    [val, firstI] = min(stepsNeeded); %the minimum number of steps needed gets there first
    
    breakdownField(NI) = breakdownField(NI) + val.*fq;
   
    if (sum(isnan(breakdownField) > 0))
        do_let_me_know = true;
    end
        % Set new BC
        
    AI = [AI; NI(firstI)];
    Occ(AI) = 1;
    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;
    
   if ~(prevAIsize == size(AI,1))
        prevAIsize = size(AI, 1);
        T = perturbation .* min_quad_with_fixed(L, [], BI, BC);
        
    end
    
    t.CData = T;
%     b.CData = breakdownField;
    s1.XData = V(AI, 1);
    s1.YData = V(AI, 2);
%     
%     s2.XData = V(AI, 1);
%     s2.YData = V(AI, 2);
    
    drawnow;
    figgif("results/DAMDiagNoisy3.gif");
    i = i + 1;
    
end


hold off