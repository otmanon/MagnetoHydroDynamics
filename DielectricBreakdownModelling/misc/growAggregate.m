function [aV, aE] = growAggregate(aV, aE, newV)
%GROWAGGREGATE Given an aggregate simplex described by aV, and aE, as well
%as new vertices, finds how to attach the new vertices to the preexisting
%aggregate. A new vertex will be connected to the closest point on the
%aggregate
numSamplesToJoin = size(newV, 1);
[D2, EI, C] = point_mesh_squared_distance(newV, aV, aE );

[aV1, aE1, nVi] = split_edges_at_point(aV, aE, EI, C);

newEEdges = [nVi', (1:numSamplesToJoin)' + size(aV1, 1)];
aV2 = [aV1; newV];
aE2 = [aE1; newEEdges];  
  
[aV3, SVI, SVJ] = remove_duplicate_vertices(aV2, 1e-5);
aE3 = SVJ(aE2);
aE4 = remove_duplicate_simplices(aE3);
aE5 = remove_degenerate_edges(aE4);
aV = aV3;
aE = aE5;
end



%Connect samples to aggregate.
%form hybrid aggregate structure
% vertInds = [1:size(aV, 1); 1:size(aV, 1)]'
% EV = [aE; vertInds];
% 
% %If closest poitn is a vertex
% Closest2Verts = EI > size(aE, 1);        %indices whrere this is 1 correspond to verts
% NC2V = sum(Closest2Verts);
% Closest2Edges = ~Closest2Verts;
% NC2E = sum(Closest2Edges);
% if (NC2V > 0)
%     newVEdges = [EI(Closest2Verts)-size(aE, 1), (1:NC2V)+ size(aV, 1)];
%     aV = [aV; randP(I(Closest2Verts), :)];
%     aE = [aE; newVEdges];
% end;

% if (NC2E > 0)
%     %If closest point is on an edge, split that edge in two. 
%     [aV, aE, nVi] = split_edges_at_point(aV, aE, EI(~Closest2Verts), C(~Closest2Verts, :));
%     newEEdges = [nVi, (1:NC2E) + size(aV, 1)];
%     aV = [aV; randP(I(~Closest2Verts), :)];
%     aE = [aE; newEEdges];  
% end;

% aE = [aE; (1:size(C, 1)) +  size(aV, 1); (1:size(randP(I, :), 1)) + size(aV, 1) + size(randP(I, :), 1) ]
% 
% aV = [aV; C; randP(I, :)];

