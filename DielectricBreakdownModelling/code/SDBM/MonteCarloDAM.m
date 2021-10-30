clc; clear;


res = 60;                  %number of grid vertices
cell_size = 1;                    %size of each cell... always make it so that the domain is the [0,1]x[0, 1] unit square

dist = 2;                   %distance from aggregate to sample from
eps = 0.1;                 %minimum radius for WoS solve.

numTimesteps = 1000;
samplingDensity = 1;
numWalks = 100;
lambda = 10;
Lb = 10;
c = 1;
numSamplesToJoin = 1;

%Initialization of grid
[nodeV, nodeF] = create_regular_grid(res);
nodeV = nodeV*res;

%Initialization of aggregate
aV =[res/2 res;res/2 res-1];        %aggregate geometry
aE = [1 2];

%Initialization of boundary
% [bV1, bE1, bN] = linspaced_square([res/2 res/2], res/2, 100);    %boundary geometry
% minY = min(bV1(:, 2));
% bottomInds = bV1(:, 2) == minY;
% bc1 = zeros(size(bV1, 1), 1);
% bc1(bottomInds) = 1000;
%
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
title('Accumulated Discharge Model');
axis([0 res 0 res]);
samplePointPlot = scatter(s1, [], [], [], [], 'filled', 'LineWidth', 0.5, 'MarkerEdgeColor', 'black');
sampledPointPlot = scatter(s1, [], [], [], 'LineWidth', 2, 'MarkerEdgeColor', 'red');

caxis([0, 1]);
newV = []; I = [];

vobj = figmp4("AccumulatedAvalanche", [])
for i=1:numTimesteps
    %update boundaryV with aggreagte
    aV1 = aV; aE1 = aE;
    bV = aV;
    bE = aE;
    
    bV = [bV1; aV];
    bE = [bE1; aE + size(bE1, 1)+1];
  
    %bc = 1000*ones(size(bV, 1), 1);                      %boundary conditions
    
    numAV = size(aV, 1);
    numBV = size(bV1, 1);
    bcA = zeros(numAV, 1);
    bc = [bc1; bcA];

    %extract level set of square distance.
    if (i == 1)
        randP = sampleVertIsosurface(aV, aE, dist, samplingDensity);
        
        phiRandP = WoS(randP, bV, bE, bc, eps, numWalks);
    
        alpha = 1/lambda * exp(- 1./(lambda*phiRandP));
        fq = c*(exp(alpha.* dist) - 1);
    else
         [randP, IJ] = sampleNewVertIsosurface(aV, aE, newV, randP, I, dist, samplingDensity);   
         newFQ = zeros(length(randP), 1);

         newFQ(1:sum(IJ)) = fq(IJ == 1);
         
         phiRandP = WoS(randP, bV, bE, bc, eps, numWalks);
             
         alpha = 1/lambda * exp(- 1./(lambda*phiRandP));
         fq = newFQ +  c*(exp(alpha.* dist) - 1);
         newV = []; I = [];
         
    end;
   
    brokenIndeces = fq >= 1.0;
    if (sum(brokenIndeces) > 0)
        I = brokenIndeces;
        newV = randP(I, :);
        [aV2, aE2] = growAggregate(aV, aE, newV);
        aV = aV2; aE = aE2;
    end;
            
    delete(aggPlotS)
    aggPlotS = plot(s1, [aV1(aE1(:, 1), 1)'; aV1(aE1(:, 2), 1)'],  [aV1(aE1(:, 1), 2)'; aV1(aE1(:, 2), 2)'], 'LineWidth', 2, 'Color', 'black');


    samplePointPlot.XData = randP(:, 1);
    samplePointPlot.YData = randP(:, 2);
    samplePointPlot.CData = fq;
    
    drawnow;
    
    if (mod(i, 5) == 0)
        vobj = figmp4("AccumulatedAvalanche", vobj);
    end;
end;



vobj.close
