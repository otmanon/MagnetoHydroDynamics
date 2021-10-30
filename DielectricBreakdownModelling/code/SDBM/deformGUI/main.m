clc; clear all; close all;

skinning_weight_types = 1;
meshPath = 'woody.obj';
[V,F] = readOBJ(meshPath);
V = V(:,1:2); % the input mesh contains redundant zero third column

% % get handles
% fig = tsurf(F,V);
% axis equal;
% fprintf( ...
%     ['Point Handle Selection: \n' ...
%     '- CLICK the mesh to add point handls \n', ...
%     '- BACKSPACE to remvoe the previous selection\n', ... 
%     '- ENTER to finish selection\n'] ...
%     );
% try
%   [Cx,Cy] = getpts;
% catch e
%   return  % quit early, stop script
% end

Cx = [156.51, 55.069, 312.21, 132.92, 242.62]';
Cy = [351.30, 243.96, 254.58, 54.055, 55.234]';
C = [Cx,Cy]; % handle locations 

% snap points to the vertices 
D = pdist2(V,C);
[~,b] = min(D); % point handle indices

% compute weights
W = compute_skinning_weight(V,F,b);


deform_GUI(V,F,C,W);
