function  FreshStableFluids()
close all; clear all;
res = 10;
h = 1/(res);
[gridV, gridF] = create_regular_grid(res);


C = cotmatrix(gridV, gridF);
G = grad(gridV, gridF);

bc = barycenter(gridV, gridF); %center of each triangle... where all vectors are stored
d = zeros(size(gridF, 1), 1);
u = zeros(size(gridF, 1), 1);
v = zeros(size(gridF, 1), 1);
div = zeros(size(gridV, 1), 1);
p = zeros(size(gridV, 1), 1);

dt = 0.1;

% [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);
% BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
hold on;
t = tsurf(gridF, gridV, 'CData', d, falpha(1, 1));
set(t, 'buttondownfcn', @onaxisdown);
colorbar();
axis equal;
colormap(parula(9));
q = [];


%do simulation here

 v = v - 0.1
while(true) 
    
    project();
    advect();
   
    t.CData = p;
    delete(q);
    q = quiver(bc(:, 1), bc(:, 2), u, v,'Color', [0, 0, 0]);
    pause(0.01);
    drawnow;
end
    function project()
        vel = reshape([u', v'], size(u, 1)*2, 1);
        div = meshDiv(gridF, vel);
        p = div\C;
        pGrad = G*p';
        pGrad = reshape(pGrad, 2, size(gridF, 1))';
        u = u - dt*pGrad(:, 1);
        v = v - dt*pGrad(:, 2);
    end

    function advect()
        vel = semiLagrangianAdvection([u, v], [u, v], dt, bc, gridV, gridF);
        u = vel(:, 1);
        v = vel(:, 2);
%         d = semiLagrangianAdvection(d, [u, v], dt, bc, gridV, gridF);
    end

    function div=meshDiv(gridF, vel);
        div = zeros(size(gridV), 1);
        div(gridF(:, 1)) = div(gridF(:, 1)) +  vel;
        div(gridF(:, 2)) = div(gridF(:, 2)) + vel;
        div(gridF(:, 3)) = div(gridF(:, 3)) + vel; 
        div = div(:, 1) + div(:, 2);
    end

    function vel = semiLagrangianAdvection(u, vel, dt, bc, V, F);
        pV = bc - dt.*vel;
        cropI = pV(:, 1) > 1 | pV(:, 1) < 0 | pV(:, 2) > 1 | pV(:, 2) < 0;
        
        [~, I, ~] = point_mesh_squared_distance(pV, V, F);
        newVel = u(I, :);
        newVel(cropI, :) = 0;     %can't get velocities from outside domain bro
        vel = newVel;
    end;
    
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