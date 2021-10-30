function AI = DAM(res, mu, sig, AI, SI);
%DAM Summary of this function goes here
%   Detailed explanation goes here
% This script implements the "Fractal Dimension of Dielectric Breakdown",
% the original DBM model introduced by N,P and W in 1984. Unphysical
% because it relies on eta parameter to make nice looking lightning bolts,
% but it is the baseline for all other models basically.

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

    i = 0;
    h = 1/res;
    lambda = 1;
    Lb = 1;
    c = 1;
    prevAIsize = 0;


    noise_mu = mu;
    noise_sig = sig;
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
    [val, firstI] = min(stepsNeeded); %the minimum number of steps needed gets there first
    
    breakdownField(NI) = breakdownField(NI) + val.*fq;
      
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
    s.XData = V(AI, 1);
    s.YData = V(AI, 2);
    
    drawnow;

    i = i + 1;
    
end


hold off
end

