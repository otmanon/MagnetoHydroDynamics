clc; clear all; close all;

meshPath = 'triwheel2D.obj';
[V,E] = readOBJ(meshPath);
%E = boundary_faces(E);
E = E(:, 1:2);
V = V(:,1:2); % the input mesh contains redundant zero third column


bI = [1; length(V)]; % first and last indices are the handles


deformLightningGUI(V,E,bI);
