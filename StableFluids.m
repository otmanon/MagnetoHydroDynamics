function  StableFluids()
close all; clear all;
res = 10;
h = 1/(res);
[gridV, gridF] = create_regular_grid(res);

bI = unique(boundary_faces(gridF));
C = cotmatrix(gridV, gridF);
G = grad(gridV, gridF);
d = zeros(size(gridV, 1), 1);
v = zeros(size(gridV, 1), 2);
g = [0, -1];
dt = 0.001;

% [~, I, ~] = point_mesh_squared_distance([rV, zeros(size(rV, 1), 1)], V3D, F);
% BaryC = barycentric_coordinates(rV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
hold on;
t = tsurf(gridF, gridV, 'CData', d, falpha(1, 1), fphong);
set(t, 'buttondownfcn', @onaxisdown);
colorbar();
axis equal;
colormap(parula(9));
q = [];


L = fdLaplacianMat(gridV, res, h);
%figure out incident faces on each vert
 

gradPV = zeros(size(gridV, 1), 2);
%do simulation here

 v = v + 0.1
while(true) 
%     v = v + g.*d;
    faceV = (1/3) * v(gridF(:, 1), :) + (1/3) * v(gridF(:, 2), :) + (1/3) * v(gridF(:, 3), :);
    faceV = reshape(permute(transpose(faceV), [2, 1]), size(gridF, 1)*2, 1);
    divV = G'*faceV;% fdDiv(v, gridV, res);

    p = -0.5*divV\L;
%     gradP = G*p';
%     gradP = transpose(permute(reshape(gradP, size(gridF, 1), 2), [2, 1]));
%     gradPV(gridF(:, 1), :) = gradP;
%     gradPV(gridF(:, 2), :) = gradP;
%     gradPV(gridF(:, 3), :) = gradP;
%     gradPV = gradPV ./3;
    gradP = fdGrad(p, gridV, res);
    v =  v - dt*gradP;
    
    
    v = semiLagrangianAdvection(v, v, dt, gridV, gridF);
    d = semiLagrangianAdvection(d, v, dt, gridV, gridF);
   
    
%      d = d + dt*L*d;
    t.CData = d;
    delete(q);
    q = quiver(gridV(:, 1), gridV(:, 2), v(:, 1), v(:, 2),'Color', [0, 0, 0]);
    pause(0.01);
    drawnow;
    max(v);
end

    function vel = semiLagrangianAdvection(u, v, dt, V, F);
        pV = V - dt.*v;
        cropI = pV(:, 1) > 1 | pV(:, 1) < 0 | pV(:, 2) > 1 | pV(:, 2) < 0;
        
        [~, I, ~] = point_mesh_squared_distance(pV, V, F);
        baryC = barycentric_coordinates(pV, V(F(I, 1), :),V(F(I, 2), :),V(F(I, 3), :));
        vel = baryC(:, 1).*u(F(I, 1), :) + baryC(:, 2).*u(F(I, 2), :) + baryC(:, 3).*u(F(I, 3), :);
        vel(cropI, :) = 0;     %can't get velocities from outside domain bro
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
    
    function gradU = fdGrad(u, gridV, res)
        u = u';
        ind = 1:size(gridV);
        [i, j] = ind2sub(res, ind);
        i = i'; j = j';
        vr = j~=res;
        vl = j~=1;
        vt = i~=res;
        vb = i~=1;
        ext = ~vr | ~vl | ~vt | ~vb;
        ri = sub2ind([res, res], i(vr), j(vr)+1);
        li = sub2ind([res, res], i(vl), j(vl)-1);
        ti = sub2ind([res, res], i(vt)+1, j(vt));
        bi = sub2ind([res, res], i(vb)-1, j(vb));
        
        
        gradU = zeros(size(gridV, 1), 2);
        
        gradU(vr, 1) = gradU(vr, 1) + u(ri);
        gradU(vl, 1) = gradU(vl, 1) - u(li);
        
        gradU(vt, 2) = gradU(vt, 2) + u(ti);
        gradU(vb, 2) = gradU(vb, 2) - u(bi);
        
        gradU(~vr, 1) = gradU(~vr, 1) + u(~vr);
        gradU(~vl, 1) = gradU(~vl, 1) - u(~vl);
        
        gradU(~vt, 1) = gradU(~vt, 2) + u(~vt);
        gradU(~vb, 1) = gradU(~vb, 2) - u(~vb);
        
        gradU = gradU * h;
        gradU(ext, :) = gradU(ext, :) * 2;
    end
    function [divU] = fdDiv(u, gridV, res)
        
        %Calculates finite difference divergence of u as specified by
        %gridded gridV
        ind = 1:size(gridV);
        [i, j] = ind2sub(res, ind);
        i = i'; j = j';
        vr = j~=res;
        vl = j~=1;
        vt = i~=res;
        vb = i~=1;
        ri = sub2ind([res, res], i(vr), j(vr)+1);
        li = sub2ind([res, res], i(vl), j(vl)-1);
        ti = sub2ind([res, res], i(vt)+1, j(vt));
        bi = sub2ind([res, res], i(vb)-1, j(vb));
        
        divU = zeros(size(gridV, 1), 1);
        divU(vr) = divU(vr) + u(ri, 1);
        divU(vl) = divU(vl) - u(li, 1);
        divU(vt) = divU(vt) + u(ti, 2);
        divU(vb) = divU(vb) - u(bi, 2);

        divU = divU * h;
        
       
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