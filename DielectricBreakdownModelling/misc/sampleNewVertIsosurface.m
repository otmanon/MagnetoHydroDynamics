function [isoP, IJ] = sampleNewVertIsosurface(V, E, newV, prevIsoP, I,  isoVal, density)
%SAMPLEISOSURFACE Samples points along the isosurface without actually
%constructing the isosurface. From randomly sampled points in your edges,
%take steps outwards along the normal of size isoVal. Then for vertices,
%sample points from the circle that is formed of radius isoVal around each
%vertex. Does one final sweep and eliminates points that are too close/far
%   INputs:
%   V   : Vertices
%   E   : Edges
%   newV: new vertices
%   prevIsoP: old points on isosurface
%   I   : Previously sampled vertex indices
%   isoVal : the value specifying the isosurface distance from the mesh
%   density : how densely do you want to sample the mesh (per unit length)


%%%%%%%%%%Now we need to handle circles around vertices

%For each vertex make a circle of radius isoVal

if (length(newV) > 0 & length(I) > 0)
    circumpherence = 2*isoVal * pi;

    numPoints = floor(circumpherence * density);   %how many points will we sample per circle
    numV = size(newV, 1);

    randVals = rand(numV*numPoints, 1)*2*pi;

    Vx = newV(:, 1);
    Vx = [repelem(Vx, numPoints, 1)];

    Vy = newV(:, 2);
    Vy = [repelem(Vy, numPoints, 1)];

    Cx = isoVal*cos(randVals);
    Cy = isoVal*sin(randVals);

    cPoints = [Vx, Vy] + [Cx, Cy];

    isoPoints2= [prevIsoP; cPoints];
    %Now find vertices that are too close and remove them
    d = point_mesh_squared_distance(isoPoints2, V, E);

    epsilon = 1e-6; %allow for floating point errors
    tooClose = (d < isoVal^2 - epsilon) | (d > isoVal^2 + epsilon);
    tooClose(I) = 1;

    IJ = ones(length(prevIsoP), 1);

    IJ(I) = 0;
    oldTooClose = (d(1:length(prevIsoP)) < isoVal^2 - epsilon) | (d(1:length(prevIsoP)) > isoVal^2 + epsilon);

    IJ(oldTooClose) = 0;

    isoP = isoPoints2(~tooClose, :);
else
    isoP = prevIsoP;
    IJ = ones(length(isoP), 1);
end;

end

