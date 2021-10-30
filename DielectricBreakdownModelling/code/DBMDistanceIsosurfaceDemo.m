% This script does the regular DBM, but instead of limiting itself to sampling grid neighbors,
% samples the distance isosurface.
clear; clc;
res = 50;
eta = 1;
[V, F] = create_regular_grid(res);
pick_max = false;
sampling_density = 1000;
iso_val = 1/res;

AV = [1 1];
AE = [1, 1];

AI = sub2ind([res, res], floor(res*AV(:, 1)), floor(res)*AV(:, 1));        %%indeces of aggregate
SI = sub2ind([res, res], floor(1), floor(1));                 %% First vertex is the source

[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI, 1), 1);
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
   % [N, NI] = getNeighbors([AX, AY], Occ);
    rV = sampleIsoSurface(AV, [], iso_val, sampling_density);
    rV = crop(rV, [0, 0], [1, 1]);
    V3D = [ V, zeros(size(V, 1), 1)];
    [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);
    
    BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
    rT = T(F(I, 1)).*BaryC(:, 1) + T(F(I, 2)).*BaryC(:, 2) + T(F(I, 3)).*BaryC(:, 3);
   % NI = snap_points(arV, V);
    

    
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
    AV = [AV; rV(RI, :)];
    [AV2, AE2] = growAggregate(AV, AE, rV(RI, :));
    AV = AV2;
    AE = AE2;
    Occ(AI) = 1;
    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;
    [BI, IA] = unique(BI);
    BC = BC(IA);
    %
    T = min_quad_with_fixed(L, [], BI, BC);
    pa = plot_edges(AV, AE, 'red', 'LineWidth', 1);
    t.CData = T;
    s.XData = V(AI, 1);
    s.YData = V(AI, 2);
    
    drawnow;
    
    
    if ( norm(rV(RI, :) - V(SI, :)) < iso_val)
        break;
    end
    
end


hold off