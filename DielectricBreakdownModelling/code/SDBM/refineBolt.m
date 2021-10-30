function [rAV, rAE] = refineBolt(V, E, rf)
%rf is the refinement factor.V, E describes bolt geometry
%   Detailed explanation goes here
assert(rf > 1);

%remove defenerate edges
degenInds = E(:, 1) == E(:, 2);
E = E(~degenInds, :);
avgD = mean(edge_lengths(V, E));
newD = avgD/rf;
eps = 1/(10*avgD);
numWalks = 100;
samplingDensity = rf*100;
rAV = [];
rAE = [];
%for each edge,refine

for e=1:size(E, 1)
    SV = V(E(e, 1), :);
    AV = V(E(e, 2), :);
    SE = [1 1];
    AE = [1 1];
    not_done = true;
    while(not_done)
        [AV, AE, newV] = MCDBMStep( AV,  AE,  SV,  SE, newD, numWalks, samplingDensity, eps);     
            
        [d, i, cp] = point_mesh_squared_distance(newV, SV,  SE);
        if (d < newD^2)
             AV = [ AV;  SV(i, :)];
             AE = [ AE; size( AV, 1)-1 size( AV, 1)];
             not_done = false;
             
             rAE = [rAE; AE + size(rAV, 1)];
             rAV = [rAV; AV];
        end
    end
end

[rAV, SVI, SVJ] = remove_duplicate_vertices(rAV, 1e-5);
rAE = SVJ(rAE);

%remove duplicate vertices


end

