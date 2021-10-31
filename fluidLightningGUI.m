
function fluidLightningGUI(V, E, sparkBI2)
close all;
N = 51;
h = 1./(N-1);

sparkV = V;
sparkE = E;
sparkBI = sparkBI2;
sparkII = setdiff(1:length(V), sparkBI);

ref_length = edge_lengths(sparkV, sparkE);
max_length = 1.5*sum(ref_length); %maximum length of entire bolt

[sparkV, sparkE] = remesh_edges(V, E, h, h/1.5);
[gridV, gridF] = create_regular_grid(N);  %first/last row/col are walls
[gridV2, gridF2] = create_regular_grid(N-2);  %first/last row/col are walls

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

div = zeros(size(gridV, 1), 1);
p = div;
default_density = 0;
d = zeros(N*N, 1);
dt = 0.01;
L = cotmatrix(gridV(bInt, :), gridF2);
hold on;
%t1 = tsurf(gridF, gridV, falpha(1, 0.0), fphong);
AEplot = [];
colorbar();
axis equal;
title("Spark Motion");
set(gca,'Color','k');
colormap(jet(101));
caxis([-1 1]);
xlim([0 1]); ylim([0 1]);



ni = [IX(floor(N/2), 4, N), IX(floor(N/2)-1, 4, N),...
    IX(floor(N/2)+1, 4, N), IX(floor(N/2)-2, 4, N),...
    IX(floor(N/2)+2, 4, N), IX(floor(N/2)-3, 4, N),...
    IX(floor(N/2)+3, 4, N), IX(floor(N/2)-4, 4, N),...
    IX(floor(N/2)+4, 4, N), IX(floor(N/2)-5, 4, N),...
    IX(floor(N/2)+5, 4, N), IX(floor(N/2)-6, 4, N),...
    IX(floor(N/2)+6, 4, N), IX(floor(N/2)-7, 4, N),...
    IX(floor(N/2)+7, 4, N)] ;

q = [];

% d(mi) =  d(mi) + 1;
% v(mi) = v(mi) + 1;
for s=1:40000
    d(IX(floor(N/2), 4, N)) = d(IX(floor(N/2), 4, N)) + 1;
    v(ni) = v(ni) + 1;
    

    project();
    handleSpark()
   
    advect();
    cellVel = interpVel(gridV(bInt, :));
    
    %
%      delete(q);
%      q = quiver(gridV(bInt, 1), gridV(bInt, 2), cellVel(:, 1), cellVel(:, 2), 1,'Color', [1, 0, 0]);

   % t1.CData = d;
    delete(AEplot);
    AEplot = plot_edges(sparkV, sparkE, 'Color', [0.9, 0.8, 0.1]);
    drawnow;
    s = s+1;
    
    if (mod(s, 2) == 0)
        figgif("spark_advection1.gif");
    end

end
    
    function handleSpark()
        advectSpark();
        [sparkV, sparkE] = remesh_edges(sparkV, sparkE, h, h/1.5);
        lengths = sum(edge_lengths(sparkV, sparkE));
        
        if (lengths > max_length) %if too long, rips and makes a new
            dist = 1/40;                   %distance from aggregate to sample from
            numWalks = 20;
            samplingDensity = 100;
            eps = 1/400;                 %minimum radius for WoS solve.
            AV = sparkV(sparkBI(1), :);
            SV = sparkV(sparkBI(2), :);
            AE = [1, 1];
            SE = [1, 1];
            not_done = true;
            while(not_done)
                [AV, AE, newV] = MCDBMStep(AV, AE, SV, SE, dist, numWalks, samplingDensity, eps);     

                [distance, i, cp] = point_mesh_squared_distance(newV, SV, SE);
                if (distance < dist^2)
                    AV = [AV; SV(i, :)];
                    AE = [AE; size(AV, 1)-1 size(AV, 1)];
                    not_done = false;
                end
            end
            bI = [1, size(AV, 1)];
            [sparkV, sparkE, sparkBI] = getDartLeader(AV, AE, bI);
            max_length = sum(edge_lengths(sparkV, sparkE))*1.5;
        
        end
    end
    function advectSpark()
        interpU = interpField(sparkV, u, [-h/2 0], N+1);
        interpV = interpField(sparkV,  v, [0 -h/2], N);
        interpD = interpField(sparkV, d, [0, 0], N);
        vel = [interpU, interpV];
        vel(sparkBI, 1:2) = 0;
        
        sparkV = sparkV + dt*vel;
        
    end

    function [x, cropI] = crop(x, maxv, minv)
        cropI = x(:, 1) >= maxv | x(:, 1) <= minv | x(:, 2) >= maxv | x(:, 2) <= minv;
        x = x(~cropI, :);
    end

    function advect()
        d = advectField(gridV, d, [0, 0], N);
        v2 = advectField(vV, v, [0 -h/2], N);
        u2 = advectField(uV, u, [-h/2 0], N+1);
        v = v2;
        u = u2;
        
    end

    function s = advectField(pos, vals, offset, NC)
        %assumes quantity d is stored at cell centers
        [pos, CI] = crop(pos, 1, 0);
        vel = interpVel(pos);
        
        prevPos = pos - dt*vel;
        [prevPos, CI2] = crop(prevPos, 1, 0);
        CI(~CI) = CI2;
        CI(CI) = 1;
        fieldVals = interpField(prevPos, vals, offset, NC);
        s = vals;
        s(~CI) = fieldVals;
        s(CI) = 0;
    end

    function outVal = interpField(pos, vals, offset, NC)
        [i, j] = world2grid(pos, h, offset);
        x = grid2world(i, j, h, offset);
        dist = pos - x;
        
        sx = dist(:, 1)/h; sy = dist(:, 2)/h;
        
        valT = (1 - sx).*vals(IX(i, j+1, NC), :) + sx.*vals(IX(i+1, j+1, NC), :);
        valB = (1 - sx).*vals(IX(i, j, NC), :) + sx.*vals(IX(i+1, j, NC), :) ;
        
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
        div = divergence(gridV);
        %poisson solve. Neumann BC on  wall vertices, divergence for rhs
        b = ones(size(L, 1), 1);
        C = dt*L;
        C = [C b; b' 0];   %bootstrap with one other equation, subtracts the constant part
        x = min_quad_with_fixed(C, [div(bInt); 0], [ni(1)], 0);
%         x = pcg(-C, [div(bInt) 0], 1e-7, 1000);
        p(bInt) = x(1:end-1);
        scale  = dt*h;
        for i = 2:N-1
            for j = 2:N-1
                v(IX(i, j+1, N)) =  v(IX(i, j+1, N)) - ...
                    scale*(p(IX(i, j+1, N)) - p(IX(i, j, N)));
                u(IX(i+1, j, N+1)) =  u(IX(i+1, j, N+1)) - ...
                    scale*(p(IX(i+1, j, N)) - p(IX(i, j, N)));
            end
        end
        [~, CIv] = crop(vV, 1-h/2, h/2);
        [~, CIu] = crop(uV, 1-h/2, h/2);
        v(CIv) = 0;
        u(CIu) = 0;
    end


    function ind = IX(i, j, NC)
        ind = i + NC*(j-1);
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


    function out = divergence(p)
        %don't do this for boundary wall cells.
        [pos, CI] = crop(p, 1, 0);
        [i, j] = world2grid(pos, h, [0, 0]);
        
        di = u(IX(i+1, j, N+1)) -  u(IX(i, j, N+1)) + ...
            v(IX(i, j+1, N)) - v(IX(i, j, N));
        
        di = -di/h;%take the negative, divide by cell spacing
        
        out = zeros(size(p, 1), 1);
        out(~CI) = di;
        
        
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