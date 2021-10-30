% This script implements the "Fractal Dimension of Dielectric Breakdown",
% the original DBM model introduced by N,P and W in 1984. Unphysical
% because it relies on eta parameter to make nice looking lightning bolts,
% but it is the baseline for all other models basically.
function [AI, AV, AE] = DBM(res, eta, AI, SI, pick_max, return_topology);

    [V, F] = create_regular_grid(res);

    [AX, AY] = ind2sub(res, AI);
    Occ = zeros(size(V, 1), 1)  ;
    Occ(AI) = 1;

    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;


    L = -cotmatrix(V, F);

    T = min_quad_with_fixed(L, [], BI, BC);

    AV = V(AI, :);
    AE = [];
    while(true)
        %Get Candidate neighbors
        [AX, AY] = ind2sub(res, AI);
        [N, NI, ANI] = getNeighbors([AX, AY], Occ);


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
        AV = [AV; V(ANI(RI), :); V(NI(RI), :)];
        AE = [AE; size(AV,1)-1 , size(AV, 1)];
        Occ(AI) = 1;
        BI = [AI; SI];
        BC = zeros(size(BI));
        BC(size(BC,1)) = 1;
        %
        T = min_quad_with_fixed(L, [], BI, BC);

    
    end
end
