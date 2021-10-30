function [V2, E2, nVI] = split_edges_at_point(V, E, I, P)
%SPLIT_EDGES_AT_POINT Splits edges of simplex given by V, E, at edge
%indeces given by I, at points aligned with I given by P.
%Inputs:
%   V - n x 2 list of vertices
%   E - m x2 list of indices indexing V
%   I - kx1 list of indices indexing E
%   P - kx2 list of points to split edges given by I around. These will be
%   new vertices
%Outpus:
%   V2  : New list of vertices
%   E2  : New list of edges
%   nVI : Indices in V2 of the newly added vertices. (May not be necessary)
%

%First things first
V2 = [V; P];
nVI= (1:size(P, 1)) + size(V, 1);

splitE = E(I, :);           %Extract edges that are undergoing splitting

En1 = splitE; En2 = splitE;
En1(:, 1) = nVI;
En2(:, 2) = nVI;



Etemp = E;
Etemp(I, :) = En1;
E2 = [Etemp; En2];



end

