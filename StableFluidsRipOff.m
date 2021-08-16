function  StableFluidsRipOff()
close all; clear all;
N = 30+2;
h = 1./(N-2);

[gridV, gridF] = create_regular_grid(N);  %first/last row/col are walls
[gridV2, gridF2] = create_regular_grid(N-2);
L = fdLaplacianMat(gridV, N, h);
bi = unique(boundary_faces(gridF));
bInt = setdiff(1:size(gridV, 1), bi);
temp = gridV(:, 1);
gridV(:, 1) = gridV(:, 2);
gridV(:, 2) = temp;
u = zeros(size(gridV, 1), 1);
v = u;
div = u;
p = u;
d = u; %ones(size(gridV, 1), 1);
dt =0.001;

C = cotmatrix(gridV, gridF);
% [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);
% BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));

% caxis([0, 10]);
s2 = subplot(3, 2, 1);
hold on;
t2 = tsurf(gridF, gridV, falpha(1, 0), fphong);
colorbar();
axis equal;
title("divergence");
colormap(parula(100));

s3 = subplot(3, 2, 2);
hold on;
t3 = tsurf(gridF, gridV, falpha(1, 0), fphong);
colorbar();
axis equal;
title("pressure");
colormap(parula(100));

s1 = subplot(3, 2, [3 4 5 6]);
hold on;
t1 = tsurf(gridF, gridV, falpha(1, 0), fphong);
colorbar();
axis equal;
title("dye");
colormap(parula(100));


% set(t, 'buttondownfcn', @onaxisdown);
mi = (N*N)/2 + N/2;
%figure out incident faces on each vert
 
% %do simulation here
 q = quiver(gridV(:, 1), gridV(:, 2), u, v, 1,'Color', [1, 1, 1]);
while(true) 
    d(mi) = d(mi) + 1;
    v(mi) = v(mi) - 1;
     project();
     advect();
    t3.CData = div;
    t2.CData = p;
    t1.CData = d;
    delete(q);
    q = quiver(gridV(:, 1), gridV(:, 2), u, v, 1,'Color', [1, 1, 1]);
    pause(0.01);
    drawnow;
    max(v);
end
    function forces();
         v(bInt) = v(bInt)  - 0.1;
         %u(bInt) = u(bInt) - 0.01.*d(bInt);
    end
    function advect();
        d = semiLagrangianAdvection(d, [u, v], dt, gridV, gridF);
        vel = semiLagrangianAdvection([u, v], [u, v], dt, gridV, gridF);
        u = vel(:, 1); v = vel(:, 2);
        
    end

    function vel = semiLagrangianAdvection(u, v, dt, V, F);
        pV = V - dt.*v;
        cropI = pV(:, 1) > 1 | pV(:, 1) < 0 | pV(:, 2) > 1 | pV(:, 2) < 0;
        
        [~, I, ~] = point_mesh_squared_distance(pV, V, F);
        baryC = barycentric_coordinates(pV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
        vel = baryC(:, 1).*u(F(I, 1), :) + baryC(:, 2).*u(F(I, 2), :) + baryC(:, 3).*u(F(I, 3), :);
        vel(cropI, :) = 0;     %can't get velocities from outside domain bro
    end;
    
    function diffuse()
        C = -cotmatrix(gridV(bInt, :), gridF2);
        M = massmatrix(gridV(bInt, :), gridF2);
        invM = 1./M;
        d = d + dt*Minv*C;
    end
    function project()
          u = set_bnd(1, u); v = set_bnd(2, v);
        %get rhs, negative divergence... but why the -0.5?
        for i = 2:N-1
            for j = 2:N-1
                div(IX(i, j)) = -0.5*h*(v(IX(i, j+1)) - v(IX(i, j-1)) ...
                    + u(IX(i+1, j)) - u(IX(i-1, j)));
                p(IX(i, j)) = 0;
            end
        end
        div = set_bnd(0, div); p = set_bnd(0, p);
        t.CData = div;
        %poisson solve. Neumann BC on  wall vertices, divergence for rhs
%          for k=1:20
%             for i = 2:N-1
%                 for j = 2:N-1
%                     p(IX(i, j)) =  (div(IX(i, j)) + p(IX(i+1, j)) + p(IX(i-1, 1))...
%                         + p(IX(i, j+1)) + p(IX(i, j-1)))/4;
%                 end
%             end
%             p = set_bnd(0, p);
%         end
        rhs = div(bInt);
        
         C = -h*h*cotmatrix(gridV(bInt, :), gridF2);
%         C = -fdLaplacianMat(gridV(bInt, :), N-2, h);
        sol = pcg(C, rhs(:), 1e-7, 3000);
        p(bInt) = sol;
        %get grad from p
        for i = 2:N-1
            for j = 2:N-1
                u(IX(i, j)) = u(IX(i, j)) - 0.5*dt*(p(IX(i+1, j)) - p(IX(i-1, j)))/(h);
                v(IX(i, j)) = v(IX(i, j)) - 0.5*dt*(p(IX(i, j+1)) - p(IX(i, j-1)))/(h);
            end
        end
        u = set_bnd(1, u); v = set_bnd(2, v);
    end
    
    function x = set_bnd(b, x)
        %set boundary conditions
        for i=2:N-1
            x(IX(1, i)) = x(IX(2, i));
            x(IX(N, i)) = x(IX(N-1, i));
            x(IX(i, 1)) = x(IX(i, 2));
            x(IX(i, N)) = x(IX(i, N-1));
            if(b==1)
                x(IX(1, i)) = -x(IX(2, i));
                x(IX(N, i)) = -x(IX(N-1, i));
            elseif(b==2)
                x(IX(i, 1)) = -x(IX(i, 2));
                x(IX(i, N)) = -x(IX(i, N-1));
            end
        end
        x(IX(1, 1)) = 0.5*(x(IX(1, 2))+ x(IX(2, 1)));
        x(IX(N, 1)) = 0.5*(x(IX(N, 2))+ x(IX(N-1, 1)));
        x(IX(1, N)) = 0.5*(x(IX(2, N))+ x(IX(1, N-1)));
        x(IX(N, N)) = 0.5*(x(IX(N, N-1))+ x(IX(N-1, N)));
    end


    function ind = IX(i, j)
        ind = i + N*(j-1);
    end;
%     function new_vel = semiLagrangianAdvection(s, vel, dt, V, F);
%         pV = V - dt.*vel;
%         cropI = pV(:, 1) > 1 | pV(:, 1) < 0 | pV(:, 2) > 1 | pV(:, 2) < 0;
%         
%         [~, I, ~] = point_mesh_squared_distance(pV, V, F);
%         baryC = barycentric_coordinates(pV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
%         new_vel = baryC(:, 1).*s(F(I, 1), :) + baryC(:, 2).*s(F(I, 2), :) + baryC(:, 3).*s(F(I, 3), :);
%         new_vel(cropI, :) = 0;     %can't get velocities from outside domain bro
%     end;
    
    function onaxisdown(src, ev);
        if (ev.Button == 1) % left click
            %  VI = [1:size(gridV,1) 1:size(gridV,1) 1:size(gridV,1)]
            %   [~, I, ~] = point_mesh_squared_distance(ev.IntersectionPoint(:, 1:2), gridV, VI );
            I = snap_points(ev.IntersectionPoint(:, 1:2), gridV);
            d(I) = 10;
        end
        
        if (ev.Button == 3) % left click
            I = snap_points(ev.IntersectionPoint(:, 1:2), gridV);
            v(I, :) = -1;
        end
    end


    function L = fdLaplacianMat(gridV, res, h)
        %0 neumann bc on the corners
           ijv = [];
           
           
           cind = (1:size(gridV, 1))';
           
           [ci, cj] = ind2sub(res, cind);
            vr = cj~=res;
            vl = cj~=1;
            vt = ci~=res;
            vb = ci~=1;
           
           ti = sub2ind([res, res],ci(vt)+1, cj(vt));
           bi = sub2ind([res, res],ci(vb)-1, cj(vb));
           ri = sub2ind([res, res],ci(vr), cj(vr)+1);
           li = sub2ind([res, res],ci(vl), cj(vl)-1);
           
           ijv = [ijv; ti, ti, -ones(length(ti), 1); ti, cind(vt), ones(length(ti), 1)];
           ijv = [ijv; bi, bi, -ones(length(bi), 1); bi, cind(vb), ones(length(bi), 1)];
           ijv = [ijv; ri, ri, -ones(length(ri), 1); ri, cind(vr), ones(length(ri), 1)];
           ijv = [ijv; li, li, -ones(length(li), 1); li, cind(vl), ones(length(li), 1)];
           
           L = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), size(gridV, 1), size(gridV, 1));
           L = L*h*h;
           
    end
    


end