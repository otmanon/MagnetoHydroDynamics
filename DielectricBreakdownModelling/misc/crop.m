function cropP = crop(P, minP, maxP)
%CROP given a list of points, gets rid of points outside of bounding box
%Inputs:
%   P:  List of points in 2D space
%   minP: bottom left coordinate of bounding box
%   maxP: top right coordinate of bounding box

croppedInds = P(:, 1) < minP(1) | P(:, 1) > maxP(1)...
           | P(:, 2) < minP(2) | P(:, 2) > maxP(2);
cropP = P(~croppedInds, :);

end

