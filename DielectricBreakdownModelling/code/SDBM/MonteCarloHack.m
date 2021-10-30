clc; clear;


res = 60;                  %number of grid vertices
cell_size = 1;                    %size of each cell... always make it so that the domain is the [0,1]x[0, 1] unit square

dist = 2;                   %distance from aggregate to sample from
eps = 0.1;                 %minimum radius for WoS solve.

numTimesteps = 300;
samplingDensity = 5;
numWalks = 5;

numSamplesToJoin = 1;

%Initialization of grid
[nodeV, nodeF] = create_regular_grid(res);
nodeV = nodeV*res;

%Initialization of aggregate
aV =[res/2 res;res/2 res-1];        %aggregate geometry
aE = [1 2];


bV1 = [res-1 0; res 0;];
bE1 = [1 2];
bc1 = [100; 100];

%%%%%%%%%%%%%%%%%%%%%% DRAWING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Position', [100, 100, 700, 600]);

s1 = subplot(1,1,1)
hold on;
aggPlotS = plot(s1, [], [], 'LineWidth', 2, 'Color', 'black');
colormap(parula(100));
colorbar();
title('Pick Max E-Field Model');
axis([0 res 0 res]);

samplePointPlot = scatter(s1, [], [], [], [], 'filled', 'LineWidth', 0.5, 'MarkerEdgeColor', 'black');
sampledPointPlot = scatter(s1, [], [], [], 'LineWidth', 2, 'MarkerEdgeColor', 'red');

vobj = figmp4("MaxEField", [])
for i=1:numTimesteps
    %update boundaryV with aggreagte
    aV1 = aV; aE1 = aE;
    
    bV = [bV1; aV];
    bE = [bE1; aE + size(bE1, 1)+1];
  
    %bc = 1000*ones(size(bV, 1), 1);                      %boundary conditions
    
    numAV = size(aV, 1);
    numBV = size(bV1, 1);
    bcA = zeros(numAV, 1);
    bc = [bc1; bcA];

    %extract level set of square distance.
    
    randP = sampleVertIsosurface(aV, aE, dist, samplingDensity);
   
    
    randP = crop(randP, [0, 0], [res, res]);

    phiRandP = WoS(randP, bV, bE, bc, eps, numWalks);

    
    
    tau = 1./ phiRandP;
    %Convert phiRandP to probability according to DBM equation
    
    [minVal, minI] = min(tau);
    minInds = tau == minVal;
    minVals = tau(minInds);
    si = randi([1 length(minVals)], 1, 1);
    
    I = find(minInds);
    I = I(si);

    
    [aV2, aE2] = growAggregate(aV, aE, randP(I, :));
    aV = aV2; aE = aE2;
    
    %Drawing info each timestep

    %drawNewAgg
   % sdfSurf.CData = nodeD;
 
    
    delete(aggPlotS)
    aggPlotS = plot(s1, [aV1(aE1(:, 1), 1)'; aV1(aE1(:, 2), 1)'],  [aV1(aE1(:, 1), 2)'; aV1(aE1(:, 2), 2)'], 'LineWidth', 2, 'Color', 'black');


    samplePointPlot.XData = randP(:, 1);
    samplePointPlot.YData = randP(:, 2);
    samplePointPlot.CData = phiRandP;
    
    drawnow;
    
    vobj = figmp4("MaxEField", vobj)
    
end;

vobj.close;



