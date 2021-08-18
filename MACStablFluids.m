function  StableFluidsRipOff()
close all; clear all;
N = 101+2;
h = 1./(N-1);

[gridV, gridF] = create_regular_grid(N);  %first/last row/col are walls

bi = unique(boundary_faces(gridF));
bInt = setdiff(1:size(gridV, 1), bi);
temp = gridV(:, 1);
gridV(:, 1) = gridV(:, 2);
gridV(:, 2) = temp;

u = zeros(N*(N+1), 1);
v = zeros(N+1 * N, 1);

div = u;
p = u;
default_density = 0;
d = u + default_density;
dt = 0.01;


L = fdLaplacianMat(gridV, N, h);

hold on;
t1 = tsurf(gridF, gridV, falpha(1, 0.0), fphong);
colorbar();
axis equal;
title("dye");
colormap(parula(100));
caxis([0 1]);
xlim([0 1]); ylim([0 1]);

% set(t, 'buttondownfcn', @onaxisdown);
mi = floor(N*N/2) + 1;
mi = [mi; mi - 1; mi + 1];

%figure out incident faces on each vert
 
% %do simulation here
%  q = quiver(gridV(:, 1), gridV(:, 2), u, v, 1,'Color', [1, 1, 1]);


% figure;
% hold on;
% t2 = tsurf(gridF, gridV, falpha(1, 0.1), fphong);
% colorbar();


for s=1:40000 
     d(mi) =  d(mi) + 0.5;
     v(mi) = v(mi) - 1;
% %       u(mi) = u(mi) - 1;
      project();
     advect();
      project();
%      delete(q);
%     q = quiver(gridV(:, 1), gridV(:, 2), u, v, 1,'Color', [1, 0, 0]);

     t1.CData = d;
%     t2.CData = p;

    drawnow;
    s = s+1;
%      if (mod(s, 166) == 0)
%          figgif("stable_fluid_longer.gif");
%      end
end
    function forces();
         v(bInt) = v(bInt)  - 0.1;
         %u(bInt) = u(bInt) - 0.01.*d(bInt);
    end
    function advect();
        [d, cropI] = semiLagrangianAdvection(d, [u, v], dt, gridV, gridF);
        d(cropI) = default_density;
        vel = semiLagrangianAdvection([u, v], [u, v], dt, gridV, gridF);
        u = vel(:, 1); v = vel(:, 2);
        
    end

    function [u2, cropI] = semiLagrangianAdvection(u1, v, dt, V, F);
        pV = V - dt.*v;
        cropI = pV(:, 1) >= 1 | pV(:, 1) <= 0 | pV(:, 2) >= 1 | pV(:, 2) <= 0;
        
        [i, j] = world2grid( pV(~cropI, :), 1/(N-1), [0, 0]);
        
        dist = pV(~cropI, :) - gridV(IX(i, j), :);
        sx = dist(:, 1)*(N-1); sy = dist(:, 2)*(N-1);
        
        valT = (1 - sx).*u1(IX(i, j+1), :) + sx.*u1(IX(i+1, j+1), :);
        valB = (1 - sx).*u1(IX(i, j), :) + sx.*u1(IX(i+1, j), :) ;
        
        val = (1 - sy).*valB + sy.*valT;
        u2 = u1;
        u2(~cropI, :) = val;
        
   end;
    
    function [i, j] = world2grid(x, h, offset)
        x = x - offset;
        ind = floor(x ./ h);
        i = ind(:, 1)+1;
        j = ind(:, 2)+1;
        
    end
        
    

    function project()
       
        %get rhs, negative divergence...
        for i = 2:N-1
            for j = 2:N-2
                div(IX(i, j)) = -0.5*(v(IX(i, j+1)) - v(IX(i, j-1)) ...
                    + u(IX(i+1, j)) - u(IX(i-1, j)))/h;
                p(IX(i, j)) = 0;
            end
        end
        
         div = set_bnd(0, div); p(:) = 0;
        %poisson solve. Neumann BC on  wall vertices, divergence for rhs

         C = -dt*L;
         p = pcg(C, div, 1e-6, 100);
%         sol = p/C;
%         for k = 1:80
%              for i = 2:N-1
%                 for j = 2:N-1
%                     p(IX(i, j)) = (div(IX(i, j)) + p(IX(i-1, j)) + p(IX(i+1, j)) +...
%                          p(IX(i, j-1)) + p(IX(i, j+1)))/4; 
%                 end
%              end
%              p = set_bnd(0, p);
%         end;
%             p = p*h*h/(dt);
        
%           p = rhs\C;
    %     p = min_quad_with_fixed(-C, div(:), [mi], [10]);

        for i = 2:N-1
            for j = 2:N-1
                u(IX(i, j)) = u(IX(i, j)) - 0.5*dt*(p(IX(i+1, j)) - p(IX(i-1, j)))/h;
                v(IX(i, j)) = v(IX(i, j)) - 0.5*dt*(p(IX(i, j+1)) - p(IX(i, j-1)))/h;
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
                x(IX(1, i)) =  -x(IX(2, i));
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
           L = L/(h*h);
           
    end
    


end