function isoP = sampleVertIsosurface(V, E, isoVal, density)
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

isoPoints2= [cPoints];
%Now find vertices that are too close and remove them
d = point_mesh_squared_distance(isoPoints2, V, E);

epsilon = 1e-6; %allow for floating point errors
tooClose = (d < isoVal^2 - epsilon) | (d > isoVal^2 + epsilon);

isoP = isoPoints2(~tooClose, :);

end

