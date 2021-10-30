function [aV, aE] = growAggregate(aV, aE, newV)
    VE = [(1:size(aV, 1))',  (1:size(aV, 1))'];
    [D2, EI, Cl] = point_mesh_squared_distance(newV, aV, VE );

    newEdge = [EI, 1 + size(aV, 1)];
    
    aV = [aV; newV];
    aE = [aE; newEdge]; 
end


