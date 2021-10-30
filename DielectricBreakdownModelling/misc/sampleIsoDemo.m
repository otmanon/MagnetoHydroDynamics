clc; clear;
%Script for sampling isosurface of a given edge.
%density
density = 10;
isoVal = 0.4;
epsilon = 0.001;

V = [1.5 1; 1.5 2; 2.5 2;];
E = [1 2; 2 3];

L = edge_lengths(V, E);
N1 = normals(V, E);
N1 = N1./vecnorm(N1, 2, 2);
N2 = -N1;

numSamplesOnCurve = floor(sum(L(:)) * 2*density);

% Get random points on curve
[randP, EInd] = random_points_on_curve(V, E, numSamplesOnCurve);

randomVals = rand(size(randP, 1), 1);
splitInds = randomVals > 0.5;

% half of the points in randP, take one normal, the other half take another
isoPoints = zeros(size(randP));
isoPoints(splitInds, :) = randP(splitInds, :) + N1(EInd(splitInds), :).*isoVal;
isoPoints(~splitInds, :) = randP(~splitInds, :) + N2(EInd(~splitInds), :).*isoVal;


%%%%%%%%%%Now we need to handle circles around vertices

%For each vertex make a circle of radius isoVal
circumpherence = 2*isoVal * pi;

numPoints = floor(circumpherence * density);   %how many points will we sample per circle
numV = size(V, 1);

randVals = rand(numV*numPoints, 1)*2*pi;

Vx = V(:, 1);
Vx = repelem(Vx, numPoints);

Vy = V(:, 2);
Vy = repelem(Vy, numPoints);

Cx = isoVal*cos(randVals);
Cy = isoVal*sin(randVals);

cPoints = [Vx, Vy] + [Cx, Cy];

isoPoints2= [isoPoints; cPoints];

d = point_mesh_squared_distance(isoPoints2, V, E);

tooClose = (d < isoVal^2 - epsilon) | (d > isoVal^2 + epsilon);

isoPoints2 = isoPoints2(~tooClose, :);
%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Position', [100 100 800 500]);

s1 = subplot(1, 2, 1);
axis([0 4 0 4]);
hold on;
plot_edges(V, E, 'LineWidth', 2);
scatter(isoPoints(:, 1), isoPoints(:, 2), 'filled', 'MarkerFaceColor', 'black');
scatter(cPoints(:, 1), cPoints(:, 2), 'filled', 'MarkerFaceColor', 'red');

s2 = subplot(1, 2, 2);
axis([0 4 0 4]);
hold on;
plot_edges(V, E, 'LineWidth', 2);
scatter(V(:, 1), V(:, 2), 'filled', 'MarkerFaceColor', 'blue');
scatter(isoPoints2(:, 1), isoPoints2(:, 2), 'filled', 'MarkerFaceColor', 'black');