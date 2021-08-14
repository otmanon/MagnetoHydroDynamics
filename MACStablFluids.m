function  FreshStableFluids()
close all; clear all; hold on;
res = 10;
h = 1/(res);

ind = 1:3;
[i, j] = meshgrid(1:res, 1:res);
i = i(:);
j = j(:);
x = g2w(i, j, [0 0]);               %i goes up in the x direction
scatter(x(:, 1), x(:, 2));          %j goes up in the y direction (THIS IS NOT row, col notation)
xlim([0 1]); ylim([0 1]);

 %center of each triangle... where all vectors are stored
[i, j] = meshgrid(1 1, 1:res);
i = i(:);
j = j(:);
x = g2w(i, j, [0 -0.5]); 
scatter(x(:, 1), x(:, 2));   
u = zeros(res, res+1);        %horizontal vels... offset by [-0.5h 0] from main grid
v = zeros(res+1, res);         %vertical vels... offset by [0 -0.5h] from main grid

u0 = zeros(res, res+1);        %horizontal vels buffer
v0 = zeros(res+1, res);         %vertical vels buffer

div = zeros(res, res);
p = zeros(res, res);
d = zeros(size(gridF, 1), 1);
dt = 0.1;

% [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);
% BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
hold on;
s = scatter();
set(t, 'buttondownfcn', @onaxisdown);
colorbar();
axis equal;
colormap(parula(9));
q = [];



%     function project()
%     end
% 
%     function advect()
%     end
% 
%     function div()
%     end

    function x = g2w(i, j, offset)
        x = [i-1, j-1] + [0.5 0.5] + offset;
        x = x*h;
    end
    function onaxisdown(src, ev);
        if (ev.Button == 1) % left click
            %  VI = [1:size(gridV,1) 1:size(gridV,1) 1:size(gridV,1)]
            %   [~, I, ~] = point_mesh_squared_distance(ev.IntersectionPoint(:, 1:2), gridV, VI );
            I = snap_points(ev.IntersectionPoint(:, 1:2), gridV);
            d(I) = 1;
            t.CData = d;
        end
        
        if (ev.Button == 3) % left click
            I = snap_points(ev.IntersectionPoint(:, 1:2), gridV);
            v(I, :) = [0, -1];
            t.CData = d;
        end
    end
    
   
end