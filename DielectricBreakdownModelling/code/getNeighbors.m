function [N, NI, AI] = getNeighbors(Coords, Occ)
    AX = Coords(:, 1);
    AY = Coords(:, 2);
    res = sqrt(size(Occ, 1));
    b = [AX - 1, AY];
    t = [AX + 1, AY];
    r = [AX, AY + 1];
    l = [AX, AY - 1];
    N = [t; r; b; l];
    ACoords = [Coords; Coords; Coords; Coords];
    AI = sub2ind([res, res], ACoords(:, 1), ACoords(:, 2));
    
    outBounds = [find(N(:, 1) > res); find(N(:, 2) > res);...
        find(N(:, 1) < 1); find(N(:, 2) < 1)];
    
    outBounds = unique(outBounds);
    inBounds = setdiff(1:size(N, 1), outBounds);
    N = N(inBounds, :);
    AI = AI(inBounds, :);
    NI = sub2ind([res, res], N(:, 1), N(:, 2));
    unOccupied = find(Occ(NI) == 0);
    
    N = N(unOccupied, :);
    NI = NI(unOccupied);
    AI = AI(unOccupied);
end
