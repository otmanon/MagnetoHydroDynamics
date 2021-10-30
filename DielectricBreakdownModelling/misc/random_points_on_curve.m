function [randP, I] = random_points_on_curve(V,E, n)
%Given 2D curve with V vertices and E edges indexing V, sample n random
%points on the curve. We want this to be a uniform sampling  along
%arclength of the curve. This amounts to sampling an edge according to a
%length weighted multinomial sampling, and then sampling a linear shape
%function within it.

%also returns I, the edge indices to which each random point falls

%Basically gptoolbox's random_points_on_mesh but for a curve.
L = edge_lengths(V, E);
L = L ./ sum(L);            %normalize edge lengths
L = [0; L];                 %WIll be endpoitns to uor binning. Need one more than we have edges

C = cumsum(L);
R = rand(n, 1);             %random points between 0, 1 used to sample edge weighted by length

[~, ~, I] = histcounts(R, C);       %Find which index R falls into. Note that this weighs edge sample by length

S = rand(n, 1);             %sample shape function

randP = [S.*V(E(I, 1), :) + (1 - S).*V(E(I, 2), :)];



end

