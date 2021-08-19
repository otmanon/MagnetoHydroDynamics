function  StableFluidsRipOff()
close all; clear all;
N = 61;
h = 1./(N-1);

[gridV, gridF] = create_regular_grid(N);  %first/last row/col are walls
temp = gridV(:, 1);
gridV(:, 1) = gridV(:, 2);
gridV(:, 2) = temp;
    

bi = unique(boundary_faces(gridF));
bInt = setdiff(1:size(gridV, 1), bi);


u = zeros(N*(N+1), 1);          % N by N+1 grid of horiz. vel 
v = zeros((N+1)*N, 1);            %N+1 by N grid of vert. vel

vV = 0;
uV = 0;
edgeGrids();                    % fills out vV and uV with the positions of the midpoint of the edges

div = u;
x = u;
default_density = 0;
d = zeros(N*N, 1);
dt = 0.01;

L = fdLaplacianMat(gridV, N, h);

hold on;
t1 = tsurf(gridF, gridV, falpha(1, 0.1), fphong);
colorbar();
axis equal;
title("dye");
colormap(parula(100));
caxis([0 1]);
xlim([0 1]); ylim([0 1]);

% set(t, 'buttondownfcn', @onaxisdown);
mi = floor(N*N/2) + 1;
mi = [mi];



%figure out incident faces on each vert
 
% %do simulation here
%  q = quiver(gridV(:, 1), gridV(:, 2), u, v, 1,'Color', [1, 1, 1]);


% figure;
% hold on;
% t2 = tsurf(gridF, gridV, falpha(1, 0.1), fphong);
% colorbar();


v = v + vV(:, 2).*vV(:, 2);
u = u + uV(:, 1).*uV(:, 1);
% u = u + 1;
q = [];
for s=1:40000 
    d(mi) =  d(mi) + 1;
%      d(mi) =  d(mi) + 0.5;
%      v(mi) = v(mi) - 1;
% %       u(mi) = u(mi) - 1;
%       project();
    
     advect();
     cellVel = interpVel(gridV(bInt, :));
     
%      
       delete(q);
      q = quiver(gridV(bInt, 1), gridV(bInt, 2), cellVel(:, 1), cellVel(:, 2), 1,'Color', [1, 0, 0]);
%         
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

   function [val, cropI] = advectVel(q, offset)
 
   end;
   
    function [x, cropI] = crop(x, maxv, minv)
        cropI = x(:, 1) >= maxv | x(:, 1) <= minv | x(:, 2) >= maxv | x(:, 2) <= minv;
        x = x(~cropI, :);
    end
    
    function advect()
        field = advectField(d);
        d(bInt) = field;
        [vx, CIx] = advectVelX();
        [vy, CIy] = advectVelY();
        u(~CIx) = vx;
        v(~CIy) = vy;
    end

    function s = advectField(d)
        %assumes quantity d is stored at cell centers
        vel = interpVel(gridV(bInt, :));
        prevPos = gridV(bInt, :) - dt*vel;
        prevPos = crop(prevPos, 1, 0);
        s = interpField(prevPos, d, [0, 0], N);
    end

    function [vx, CI]=advectVelX()
        %For each vertical edge  (Nx(N+1))
        
        %must interpolate vy
        [pos, CI] = crop(uV, 1, 0);
        vel = interpVel(pos);
        
        prevPos = pos - dt*vel;
        [prevPos, CI2] = crop(prevPos, 1, 0);
        CI(~CI) = CI2;
        CI(CI) = 1;
        vx = interpField(prevPos, u, [-h/2 0], N+1);
    end

    function [vy, CI] = advectVelY()
       
        %must interpolate vy
        [pos, CI] = crop(vV, 1, 0);
        vel = interpVel(pos);
        
        prevPos = pos - dt*vel;
        [prevPos, CI2] = crop(prevPos, 1, 0);
        CI(~CI) = CI2;
        CI(CI) = 1;
        vy = interpField(prevPos, v, [0 -h/2], N);
        
    end

    function outVal = interpField(p, val, offset, NC)
        [i, j] = world2grid(p, h, offset);
        x = grid2world(i, j, h, offset);
        dist = p - x;
        
        sx = dist(:, 1)/h; sy = dist(:, 2)/h;
        
        valT = (1 - sx).*val(IX(i, j+1, NC), :) + sx.*val(IX(i+1, j+1, NC), :);
        valB = (1 - sx).*val(IX(i, j, NC), :) + sx.*val(IX(i+1, j, NC), :) ;
        
        outVal = (1 - sy).*valB + sy.*valT;
    end
    function vel=interpVel(p)
        velY = interpField(p, v, [0 -h/2], N ); 
        velX = interpField(p, u, [-h/2 0], N+1 );
        vel = [velX velY];
    end
   

    function p = grid2world(i, j, h, offset)
        p = [i j] - 1;
        p = p .* h;
        p = p + offset;
    end

    function [i, j] = world2grid(x, h, offset)
        x = x - offset;
        ind = floor(x ./ h);
        i = ind(:, 1) + 1;
        j = ind(:, 2) + 1;
    end
        
    

    function project()
       
        %get rhs, negative divergence...
        for i = 2:N-1
            for j = 2:N-2
                div(IX(i, j)) = -0.5*(v(IX(i, j+1)) - v(IX(i, j-1)) ...
                    + u(IX(i+1, j)) - u(IX(i-1, j)))/h;
                x(IX(i, j)) = 0;
            end
        end
        
         div = set_bnd(0, div); x(:) = 0;
        %poisson solve. Neumann BC on  wall vertices, divergence for rhs

         C = -dt*L;
         x = pcg(C, div, 1e-6, 100);
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
                u(IX(i, j)) = u(IX(i, j)) - 0.5*dt*(x(IX(i+1, j)) - x(IX(i-1, j)))/h;
                v(IX(i, j)) = v(IX(i, j)) - 0.5*dt*(x(IX(i, j+1)) - x(IX(i, j-1)))/h;
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


    function ind = IX(i, j, NC)
        ind = i + NC*(j-1);
    end;
    
    function ind = IXv(i, j)
        ind = i + (N+1)*(j-1);
    end;
    
    function ind = IXu(i, j)
        ind = i + (N)*(j-1);
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
    

    function edgeGrids()
        %get grid positions of edge midpoints. Assume cell centers go from
        %(0, 0) to (1, 1). 
        uV = zeros(size(u, 1), 2);
        ind = 1;
        %build u grid
        for j=1:N
            for i=1:N+1
                uV(ind, :) = grid2world(i, j, h, [-h/2 0]); 
                ind = ind + 1;
            end
        end
        
        ind = 1;
        vV = zeros(size(v, 1), 2);
        for j=1:N+1
            for i=1:N
                vV(ind, :) = grid2world(i, j, h, [0 -h/2]);
                ind = ind + 1;
            end
        end
        %build v grid
    end;


end