function [dartV, dartE, dartBI] = getDartLeader(V, E, bI)
%extracts shortest path to bI 
    assert(length(bI) == 2);
    A = adjacency_matrix(E);
    G = graph(A);
    
    path = shortestpath(G, bI(1), bI(2));
    
    %dartV = V(path, :);
    dartE = [path(1:end-1)' path(2:end)'];
    
    [dartV, I] = remove_unreferenced(V, dartE);
    dartE = I(dartE);
    dartBI = I(bI);
end

