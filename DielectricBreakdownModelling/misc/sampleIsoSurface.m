function isoP = sampleIsoSurface(V, E, isoVal, density)
%SAMPLEISOSURFACE Samples points along the isosurface without actually
%constructing the isosurface. From randomly sampled points in your edges,
%take steps outwards along the normal of size isoVal. Then for vertices,
%sample points from the circle that is formed of radius isoVal around each
%vertex. Does one final sweep and eliminates points that are too close/far
%   INputs:
%   V   : Vertices
%   E   : Edges
%   isoVal : the value specifying the isosurface distance from the mesh
%   density : how densely do you want to sample the mesh (per unit length)

isoPoints = [];
if (size(E, 1) > 0)
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
else
    E = [[1:size(V, 1)]', [1:size(V, 1)]'];
end;
 

%%%%%%%%%%Now we need to handle circles around vertices

%For each vertex make a circle of radius isoVal
circumpherence = 2*isoVal * pi;

numPoints = floor(circumpherence * density);   %how many points will we sample per circle
numV = size(V, 1);

randVals = rand(numV*numPoints, 1)*2*pi;

Vx = V(:, 1);
Vx = [repelem(Vx, numPoints, 1)];

Vy = V(:, 2);
Vy = [repelem(Vy, numPoints, 1)];

Cx = isoVal*cos(randVals);
Cy = isoVal*sin(randVals);

cPoints = [Vx, Vy] + [Cx, Cy];

isoPoints2= [isoPoints; cPoints];

%Now find vertices that are too close and remove them
d = point_mesh_squared_distance(isoPoints2, V, E);

epsilon = 1e-6; %allow for floating point errors
tooClose = (d < isoVal^2 - epsilon) | (d > isoVal^2 + epsilon);

isoP = isoPoints2(~tooClose, :);

end

