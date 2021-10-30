% This script , instead of sampling a distance isosurface, samples the
% isosurface of the laplace field.
clear; clc;
res = 100;
eta = 1;
[V, F] = create_regular_grid(res);
pick_max = false;
sampling_density = 100;
iso_val = 8/res;
AV = [0.75, 0.75];
AE = [1, 1];
AI = sub2ind([res, res], floor(res*AV(1)), floor(res*AV(2)));        %%indeces of aggregate
SI = sub2ind([res, res], floor(res/2), floor(res/2));                 %% First vertex is the source


[AX, AY] = ind2sub(res, AI);
Occ = zeros(size(V, 1), 1)  ;
Occ(AI) = 1;

BI = [AI; SI];
BC = zeros(size(BI));
BC(size(BC,1)) = 1;


L = -cotmatrix(V, F);

%implement periodic Bc
topi = find(V(:, 2) == max(V(:, 2)));
boti = find(V(:, 2) == min(V(:, 2)));
lefti =  find(V(:, 1) == min(V(:, 1)));
righti = find(V(:, 1) == max(V(:, 1)));
numB = size(topi, 1) + size(righti, 1);

Aeq = sparse(numB, size(V, 1));
topDiag = sub2ind([numB, size(V, 1)], [1:size(topi, 1)]', topi);
rightDiag = sub2ind([numB, size(V, 1)], [size(topi, 1)+1:numB]', righti);
topOffDiag = sub2ind([numB, size(V, 1)], [1:size(topi, 1)]', boti);
rightOffDiag = sub2ind([numB, size(V, 1)], [size(topi, 1)+1:numB]', lefti);
Aeq(topDiag) = 1;
Aeq(rightDiag) = 1;
Aeq(topOffDiag) = -1;
Aeq(rightOffDiag) = -1;
T = min_quad_with_fixed(L, [], BI, BC, Aeq, []);


%get laplacian isolines
[LS, LD, I] = isolines(V, F, T, iso_val);
isoV = [LS;LD];
isoE = [ [1:size(LS, 1)]', size(LS, 1) + [1:size(LD, 1)]'];
[isoV, SVI, SVJ] = remove_duplicate_vertices(isoV, 1e-6);
isoE = SVJ(isoE);


fig = figure('Position', [100, 100, 1200, 1000]);
hold on;
t = tsurf(F, V, 'CData', T, fphong, falpha(1, 0));
s = scatter(V(AI, 1), V(AI, 2), 'filled', 'MarkerFaceColor', 'red');
p = plot_edges(isoV, isoE, 'black');
%pa = plot([LS(:, 1)'; LD(:, 1)'],[LS(:, 2)'; LD(:, 2)'], 'black', 'LineWidth', 1);

colorbar
while(true)
    %Get Candidate neighbors
    [AX, AY] = ind2sub(res, AI);
    
    %get laplacian isoline
    [LS, LD, I] = isolines(V, F, T, iso_val);
    isoV = [LS;LD];
    isoE = [ [1:size(LS, 1)]', size(LS, 1) + [1:size(LD, 1)]'];
    [isoV, SVI, SVJ] = remove_duplicate_vertices(isoV, 1e-6);
    isoE = SVJ(isoE);
    
    rV = random_points_on_curve(isoV, isoE, 1);
   % [AV2, AE2] = growAggregate(AV, AE, rV);
   % AV = AV2; AE = AE2;
    NI = snap_points(rV, V);
    % Set new BC
    AI = [AI; NI];
    BI = [AI; SI];
    BC = zeros(size(BI));
    BC(size(BC,1)) = 1;
    [BI, IA] = unique(BI);
    BC = BC(IA);
    %
    T = min_quad_with_fixed(L, [], BI, BC, Aeq, []);
    
    t.CData = T;
    s.XData = V(AI, 1);
    s.YData = V(AI, 2);
    p = plot_edges(isoV, isoE, 'black');
   % pa = plot_edges(AV, AE, 'red', 'LineWidth', 2);
    drawnow;
    
    figgif("./results/samplingLaplacianIsocontourPeriodicBC.gif");
    
    if ( norm(rV - V(SI, :)) < iso_val)
        break;
    end
    
end


hold off