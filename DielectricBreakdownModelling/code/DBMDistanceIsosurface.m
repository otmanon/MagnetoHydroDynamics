% This script does the regular DBM, but instead of limiting itself to sampling grid neighbors,
% samples the distance isosurface.
function [AI, AV] = DBMDistanceIsoSurface(res, eta, AI, SI, pick_max, sampling_density)

    [V, F] = create_regular_grid(res);

    iso_val = 1/res;

    [AX, AY] = ind2sub(res, AI);
    Occ = zeros(size(V, 1), 1)  ;
    Occ(AI) = 1;

    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;


    L = -cotmatrix(V, F);

    T = min_quad_with_fixed(L, [], BI, BC);

    while(true)
        %Get Candidate neighbors
        [AX, AY] = ind2sub(res, AI);
       % [N, NI] = getNeighbors([AX, AY], Occ);
        rV = sampleIsoSurface(V(AI, :), [], iso_val, sampling_density);
        rV = crop(rV, [0, 0], [1, 1]);
        V3D = [ V, zeros(size(V, 1), 1)];
        [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);

        BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
        rT = T(F(I, 1)).*BaryC(:, 1) + T(F(I, 2)).*BaryC(:, 2) + T(F(I, 3)).*BaryC(:, 3);
       % NI = snap_points(rV, V);

    
        if (sum(rT) == 0)
            break;
        end
        W = rT;
        W = W .^ eta;
        W = W ./ sum(W);
        %Sample one of them
        RI = weightedSampling(W, 1);
        if pick_max == true
            [~, RI] = max(W);
        end

        
        NI = snap_points(rV(RI, :), V);
        % Set new BC
        AI = [AI; NI];
        Occ(AI) = 1;
        BI = [AI; SI];
        BC = zeros(size(BI));
        BC(size(BC,1)) = 1;
        [BI, IA] = unique(BI);
        BC = BC(IA);
        %
        T = min_quad_with_fixed(L, [], BI, BC);




        if ( norm(rV(RI, :) - V(SI, :)) < iso_val)
            break;
        end

    end
    
    AV = V(AI, :);

end