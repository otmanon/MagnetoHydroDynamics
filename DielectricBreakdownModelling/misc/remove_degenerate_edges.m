function E2 = remove_degenerate_edges(E)
%goes through edge list. If an edge provides a self-loop (both vertex
%indices are the same) , remove it.
    ind1 = E(:, 1);
    ind2 = E(:, 2);
    degInd = ind1 == ind2;
    E2 = E(~degInd, :);
end
