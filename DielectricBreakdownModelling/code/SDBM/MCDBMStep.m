function [AV, AE, newV] = MCDBMStep(AV, AE, SV, SE, dist, numWalks, samplingDensity, eps )
%MCDBM Summary of this function goes here
%   Detailed explanation goes here
    bV = [ AV; SV];
    bE = [ AE; SE+size(AV, 1)];
    bc = [zeros(size( AV, 1), 1);...
        ones(size( SE, 1), 1)];              %zeros at aggregate, ones at source

    randP = sampleVertIsosurface( AV,  AE, dist, samplingDensity);

    phiRandP = WoS(randP, bV, bE, bc, eps, numWalks);

    tau = 1./ phiRandP;
    [minVal, minI] = min(tau);
    minInds = tau == minVal;
    minVals = tau(minInds);
    si = randi([1 length(minVals)], 1, 1);
    I = find(minInds);
    I = I(si);
    newV = randP(I, :);
    [AV2, AE2] = growAggregateVert( AV,  AE, newV);
     AV = AV2;  AE = AE2;
end

