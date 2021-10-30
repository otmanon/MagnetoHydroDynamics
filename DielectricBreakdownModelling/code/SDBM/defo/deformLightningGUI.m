function deformLightningGUI(V,E,bI)
mainfig = figure
fprintf( ...
    ['Linear Blend Skinning: \n' ...
    '- CLICK handles to visualize weights \n', ...
    '- DRAG a handle to move\n', ...
    '- SHIFT+DRAG a handle to rotate\n'] ...
    );
C = V(bI, :);
DeformOBJ.new_C = C;
np = numel(bI);  % number of point handles
P = 1:size(C,1);
DeformOBJ.bI = bI;
% keep track of rotations stored at each control point, for 2D this is a m
% by 1 list of angles
DeformOBJ.R = zeros(np,1);
DeformOBJ.update_positions = @update_positions;


%DeformOBJ.tsh = plot_edges(V, E, 'Color', [1 0 0], );;
axis equal
hold on
set(gca,'ButtonDownFcn',@onaxisdown)
C_plot = scatter3( ...
    C(:,1),C(:,2),0.1+0*C(:,1), ...
    'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
    'LineWidth',2,'SizeData',100, ...
    'ButtonDownFcn',@oncontrolsdown);
set(gca,'Color','k')

% set(fig,'KeyPressFcn',@onkeypressdown);

% axis manual


% plot the original mesh
DeformOBJ.wvsh = plot_edges(V, E, 'Color', [1 1 0]);
DeformOBJ.VD = V;
DeformOBJ.V = V;
DeformOBJ.E = E;
view(2);
xlim([0 1]); ylim([0 1]);

% axis manual

%Get all the differential operators like stiffness/laplacian matrix



%Curve Laplacian Matrix
L = curveLaplacian(DeformOBJ.V, DeformOBJ.E);
K = assembleK(DeformOBJ.V, DeformOBJ.E);
lengths = edge_lengths(DeformOBJ.V, DeformOBJ.E);



title('Deformable Lightning');


% keep track of window xmin, xmax, ymin, ymax



% keep track of down position

alpha=0;
down_pos = [];
% keep track of last two drag positions
drag_pos = [];
last_drag_pos = [];
% keep track of mesh vertices at mouse down
down_V = [];
% keep track of index of selected control point
ci = [];
% type of click ('left','right')
down_type  = '';

EHistDelta = inf;

% Callback for setting new control points on mesh
    function onaxisdown(src, ev)
        if(ev.Button == 1)
            % left-click
            down_type = 'left';
        elseif(ev.Button == 2)
            % other (right) click
            down_type = 'right';
        end
        down_pos=get(gca,'currentpoint');
        [minD,vi] =  ...
            min(sum((DeformOBJ.VD(:,1:2) - ...
            repmat(down_pos(1, 1:2),size(DeformOBJ.VD,1),1)).^2,2));
        [found, iB] = ismember(vi, DeformOBJ.bI);
        
        if (~found(1) && strcmp(down_type, 'left'))
            DeformOBJ.bI = [DeformOBJ.bI; vi];
            DeformOBJ.new_C = [DeformOBJ.new_C; DeformOBJ.VD(vi, :)];
            
            delete(C_plot)
            C_plot = scatter3( ...
                DeformOBJ.new_C(:,1),DeformOBJ.new_C(:,2),0.1+0*DeformOBJ.new_C(:,1), ...
                'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
                'LineWidth',2,'SizeData',100, ...
                'ButtonDownFcn',@oncontrolsdown);
        end;
        
    end

% Callback for mouse down on control points
    function oncontrolsdown(src,ev)
        
        drag_pos=get(gca,'currentpoint');
        drag_pos=[drag_pos(1,1,1) drag_pos(1,2,1)];
        % get current mouse position, and remember old one
        down_pos=get(gca,'currentpoint');
        down_pos=[down_pos(1,1,1) down_pos(1,2,1)];
        last_drag_pos=down_pos;
        drag_pos=down_pos;
        % keep track of control point positions at mouse down
        DeformOBJ.new_C = [get(C_plot,'XData')' get(C_plot,'YData')'];
        % get index of closest control point
        [minD,ci] =  ...
            min(sum((DeformOBJ.new_C(:,1:2) - ...
            repmat(down_pos(:, 1:2),size(DeformOBJ.new_C,1),1)).^2,2));
        % keep track of mesh vertices at mouse down
        down_V = DeformOBJ.V;
        
        % tell window that drag and up events should be handled by controls
        set(gcf,'windowbuttonmotionfcn',@oncontrolsdrag)
        set(gcf,'windowbuttonupfcn',@oncontrolsup)
        set(gcf,'KeyPressFcn',@onkeypress)
        if(ev.Button == 1)
            % left-click
            down_type = 'left';
        elseif(ev.Button == 3)
            % other (right) click
            down_type = 'right';
        end
        
        [minD,vi] =  ...
            min(sum((DeformOBJ.VD(:,1:2) - ...
            repmat(down_pos(1, 1:2),size(DeformOBJ.VD,1),1)).^2,2));
        [found, iB] = ismember(vi, DeformOBJ.bI);
        if (found(1) && strcmp(down_type, 'right'))
            ind = 1:length(DeformOBJ.bI);
            DeformOBJ.bI = DeformOBJ.bI(ind ~= iB); % delete bI
            DeformOBJ.new_C = DeformOBJ.new_C(ind ~= iB, :);
            
            delete(C_plot)
            C_plot = scatter3( ...
                DeformOBJ.new_C(:,1),DeformOBJ.new_C(:,2),0.1+0*DeformOBJ.new_C(:,1), ...
                'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
                'LineWidth',2,'SizeData',100, ...
                'ButtonDownFcn',@oncontrolsdown);
            update_positions();
            %get new positions for all other verts
            
            deform_curve_ARAP();
            update_deformed_positions();
        end;
        % try to find ci in list of point handles
        [found, iP] = ismember(ci,P);
        
    end

% Callback for mouse drag on control points
    function oncontrolsdrag(src,ev)
        % keep last drag position
        % get current mouse position
        drag_pos=get(gca,'currentpoint');
        drag_pos=[drag_pos(1,1,1) drag_pos(1,2,1)];
        
        last_drag_pos = DeformOBJ.new_C(ci, 1:2);
        drag = drag_pos - last_drag_pos;
        if(strcmp('left',down_type))
            % move selected control point by drag offset
            DeformOBJ.new_C(ci,:) = ...
                DeformOBJ.new_C(ci,1:2) + drag_pos(:, 1:2)-last_drag_pos(:, 1:2);
            
        else
            [found, iP] = ismember(ci,P);
            if(found)
                DeformOBJ.R(iP) = ...
                    DeformOBJ.R(iP) + 2*pi*(drag_pos(1)-last_drag_pos(1))/100;
            end
        end
        update_positions();
        
        stretch = deform_curve_ARAP();
        %get new positions for all other verts
        while (max(stretch) > 1000000)
            dielectric_breakdown_step();
            L = curveLaplacian(DeformOBJ.V, DeformOBJ.E);
            K = assembleK(DeformOBJ.V, DeformOBJ.E);
            lengths = edge_lengths(DeformOBJ.V, DeformOBJ.E);
            % delete(DeformOBJ.tsh)
            % DeformOBJ.tsh = plot_edges(DeformOBJ.V, DeformOBJ.E, 'Color', [1 0 0], );;
            stretch = deform_curve_ARAP();
        end;
        
        update_deformed_positions();
        
        
    end

    function update_deformed_positions()
        
        delete(DeformOBJ.wvsh)
        DeformOBJ.wvsh = plot_edges(DeformOBJ.VD, DeformOBJ.E, 'Color', [1 1 0] );
        
    end;
    
    
    function stretch = deform_curve_ARAP()
        
        bI = DeformOBJ.bI;
        bV = DeformOBJ.new_C;
        V0 = DeformOBJ.V;
        
        %      bending_energy_quad = alpha*L'*invM*L;
        %      bending_energy_linear = -2*alpha*L'*invM*L*V0;
        % %
        %      stretching_energy_quad = beta*0.5*k*G;
        
        bending_energy_quad = -L;
        precompute = [];
        x = DeformOBJ.V;
               
        for i=1:20
            Covs = x'*K;
            Clist = reshape(Covs, 2, 2, size(DeformOBJ.V, 1));
            Rlist = fit_rotations(permute(Clist, [2, 1, 3]));
            Rstack = transpose(reshape(Rlist, 2, 2*size(DeformOBJ.V, 1)));
            bending_energy_linear = -K*Rstack;
            [x, precompute] = min_quad_with_fixed(bending_energy_quad, ...
                bending_energy_linear, bI, bV);
        end;
        
        DeformOBJ.VD = x;
        stretch = edge_lengths(DeformOBJ.VD, DeformOBJ.E)./edge_lengths(DeformOBJ.V,  DeformOBJ.E);
        
        
        
        
        
    end;
    
    function dielectric_breakdown_step()
        
        bV1 = [DeformOBJ.new_C(ci, :); DeformOBJ.new_C(ci, :) + 0.01];
        bE1 = [(1:size(DeformOBJ.new_C(ci, :),1))', size(DeformOBJ.new_C(ci, :), 1) + (1:size(DeformOBJ.new_C(ci, :),1))' ];
        bc1 = ones(2*size(DeformOBJ.new_C(ci, :), 1), 1);
        aV = DeformOBJ.V;
        aE = DeformOBJ.E;
        
        
        %update boundaryV with aggreagte
        
        bV = [bV1; aV];
        bE = [bE1; aE + size(bV1, 1)];
        
        numAV = size(aV, 1);
        numBV = size(bV1, 1);
        bcA = zeros(numAV, 1);
        bc = [bc1; bcA];
        
        %extract level set of square distance.
        dist = mean(edge_lengths(DeformOBJ.V, DeformOBJ.E));
        samplingDensity=100;
        randP = sampleVertIsosurface(aV, aE, dist, samplingDensity);
        
        numWalks = 5;
        phiRandP = WoS(randP, bV, bE, bc, 1e-7, numWalks);
        
        
        
        tau = 1./ phiRandP;
        %Convert phiRandP to probability according to DBM equation
        
        [minVal, minI] = min(tau);
        minInds = tau == minVal;
        minVals = tau(minInds);
        si = randi([1 length(minVals)], 1, 1);
        
        I = find(minInds);
        I = I(si);
        
        
        [aV2, aE2] = attachNewStrike(aV, aE, randP(I, :));
        aV = aV2; aE = aE2;
        DeformOBJ.V = aV;
        DeformOBJ.E = aE;
    end

    function [aV, aE] = attachNewStrike(aV, aE, newV)
        
        VE = [(1:size(aV, 1))',  (1:size(aV, 1))'];
        [D2, EI, Cl] = point_mesh_squared_distance(newV, aV, VE );
        
        newEdge = [EI, 1 + size(aV, 1)];
        [found, index] = ismember(EI, DeformOBJ.bI);
        if (found(1))
            DeformOBJ.bI(index) = 1 + size(aV, 1);
        end
        
        aV2 = [aV; newV];
        aE2 = [aE; newEdge];
        
        aV = aV2;
        aE = aE2;
    end
    function update_positions()
        % update display positions
        
        delete(C_plot)
        C_plot = scatter3( ...
            DeformOBJ.new_C(:,1),DeformOBJ.new_C(:,2),0.1+0*DeformOBJ.new_C(:,1), ...
            'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
            'LineWidth',2,'SizeData',100, ...
            'ButtonDownFcn',@oncontrolsdown);
        
        % USING LINEAR BLEND SKINNING
        % get transformations stored at each point and bone handle
        
        
        % update mesh positions
    end

% Callback for mouse release of control points
    function oncontrolsup(src,ev)
        % Tell window to handle drag and up events itself
        set(gcf,'windowbuttonmotionfcn','');
        set(gcf,'windowbuttonupfcn','');
        cur_V = DeformOBJ.V;
        
        % scale window to fit
        
    end

    function onkeypress(src,ev)
        if(strcmp(ev.Character,'r'))
            DeformOBJ.bI = bI;
            DeformOBJ.new_C = DeformOBJ.V(DeformOBJ.bI, :);
            DeformOBJ.VD = DeformOBJ.V;
            update_positions();
            update_deformed_positions();
            
        elseif(strcmp(ev.Key,'uparrow'))
            
            alpha = min(alpha + 0.1, 1.0);
            beta = 1 - alpha;
            title(strcat('α : ', num2str(alpha)));
            
        elseif(strcmp(ev.Key, 'downarrow'))
            alpha = max(alpha - 0.1, 0.0);
            beta = 1 - alpha;
            title(strcat('α : ', num2str(alpha)));
        elseif(strcmp(ev.Character, 'a'))
            %animate with circular motion at control points
            animateControlPoints();
        elseif(strcmp(ev.Character, 'f'))
            [dartV, dartE, dartBI] = getDartLeader(DeformOBJ.VD, DeformOBJ.E, DeformOBJ.bI);
            
            fluidLightningGUI(dartV, dartE);
        end
        
    end


    
    function animateControlPoints()
        %get circle for each control point
        
        animfig = figure;
        axis equal
        hold on
        xlim([0 1]); ylim([0, 1]);
        set(gca,'Color','k')
        title('Lightning Animation');
        view(2);
        %rad of each circle
        r = 0.01*ones(size(DeformOBJ.new_C, 1), 1) + 0.01*rand(size(DeformOBJ.new_C, 1), 1);
        
        %get a random angle between 0 and 2pi
        theta = 2*pi*rand(size(DeformOBJ.new_C, 1), 1);
        speed = 2*rand(size(DeformOBJ.new_C, 1), 1);
        ccw = 2*(round(rand(size(DeformOBJ.new_C, 1), 1)) - 0.5);
        anim_res = 100;
        theta_step = 2*pi/anim_res;
        
        dt = 0.01;
       
        figure(animfig);
        

        
        [dartV, dartE, dartBI] = getDartLeader(DeformOBJ.VD, DeformOBJ.E, DeformOBJ.bI);
        Mmat = massmatrix(dartV, dartE);
        Lmat = curveLaplacian(dartV, dartE);
        Kmat = assembleK(dartV, dartE);
        c = 0.01;
        k = 0.01;
        qn = dartV;
        qnm1 = dartV;
        qnp1 = []; %solving for this.
        for i=5:anim_res
            y = 2*qn - qnm1;
           % y = qn; with non inertial forces
            quad_energy = 0.5*Mmat - dt*dt*c*Lmat - dt*dt*k*Lmat;
            linear_energy = -Mmat*y;
            
            x = dartV;
            for i=1:20
                Covs = x'*Kmat;
                Clist = reshape(Covs, 2, 2, size(dartV, 1));
                Rlist = fit_rotations(permute(Clist, [2, 1, 3]));
                Rstack = transpose(reshape(Rlist, 2, 2*size(dartV, 1)));
                linear_energy_real = linear_energy-dt*dt*k*Kmat*Rstack;
                [x, precompute] = min_quad_with_fixed(quad_energy, ...
                    linear_energy_real, dartBI, DeformOBJ.new_C);
            end;
        
            qnp1 = x;
            qnm1 = qn;
            qn = qnp1;
            delete(DeformOBJ.wvsh)
            DeformOBJ.wvsh = plot_edges(qnp1, dartE, 'Color', [1 1 0] );
            drawnow;
          
            %           a
            vx = speed.*r.*cos(theta);
            vy = speed.*r.*sin(ccw.*theta);
            vel = [vx vy];
            DeformOBJ.new_C(2, :) = DeformOBJ.new_C(2, :) + 0.2*vel(2, :);
            
            
            theta = theta + theta_step;
            if (mod(i, 2) == 0)
                figgif('results/mcfARAPinertialLightningC0.01K0.01.gif');
            end;
        end
         figure(mainfig)
         update_positions();
         update_deformed_positions();
         drawnow();
        
    end

    

end