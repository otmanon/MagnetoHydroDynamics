function lightningGrowthUI()
clc; clear all; close all;
fprintf( ...
    ['Lightning Growth App: \n' ...
    '- LEFT CLICK to set a starting point \n', ...
    '- RIGHT CLICK to set a source/destination\n', ...
    ] ...
    );

App.init = 'true';
App.growthAlg = 'IGDBM';
App.growthAlg = 'MCDBM';
App.mouse_up = false;
App.not_done = true;
App.AEplot = [];
App.AV = []; App.AE = []; App.SV = []; App.SE = [];
App.AV0 = []; App.AE0 = []; App.SV0 = []; App.SE0 = [];

hold on;
axis equal;
set(gca,'ButtonDownFcn',@mouse_click_callback)
set(gca,'Color','k')
set(gcf,'KeyPressFcn',@keyboard_down_callback)
xlim([0, 1]); ylim([0, 1]);
title("Stepped Leader Growth");

    function mouse_click_callback(src, ev )
        App.down_pos=get(gca,'currentpoint');
        if(ev.Button == 1)
            % left-click
            App.down_type = 'left';
            App.AV = [App.AV; App.down_pos(1, 1:2)]; App.AV0 = [App.AV0; App.down_pos(1, 1:2)] ;
            App.AE = [App.AE; size(App.AV, 1) size(App.AV, 1)]; App.AE0 = [App.AE0; size(App.AV0, 1) size(App.AV0, 1)];
            App.controlAVplot = scatter(App.AV(:, 1), App.AV(:, 2), 'filled',...
            'MarkerFaceColor',[0.9 0.8 0.1]);
        elseif(ev.Button == 3)
            % right click
            App.down_type = 'right';
            App.SV = [App.SV; App.down_pos(1, 1:2)]; App.SV0 = [App.SV0; App.down_pos(1, 1:2)] ;
            App.SE = [App.SE; size(App.SV, 1) size(App.SV, 1)]; App.SE0 = [App.SE0; size(App.SV0, 1) size(App.SV0, 1)];
            App.controlSVplot = scatter(App.SV(:, 1), App.SV(:, 2),  'filled',...
            'MarkerFaceColor',[0.9 0.1 0.1]);
        end 
 

    end

    function keyboard_down_callback(src, ev)
        if(strcmp(ev.Character,'r'))
            App.AV = App.AV0;
            App.AE = App.AE0;
            App.SV = App.SV0;
            App.SE = App.SE0;
            delete(App.AEplot);
        elseif(strcmp(ev.Character, 'c'))
            App.AV = [];
            App.AE = [];
            App.SV = [];
            App.SE = [];
            App.AV0 = [];
            App.AE0 = [];
            App.SV0 = [];
            App.SE0 = [];
            delete(App.controlSVplot);
            delete(App.controlAVplot);
            delete(App.AEplot);
        elseif(strcmp(ev.Key, '1'))
            fprintf(['Monte Carlo Algorithm']);
            title("Stepped Leader Growth: Monte Carlo DBM");
            App.growthAlg = 'MCDBM';
        elseif(strcmp(ev.Key, '2'))
            fprintf(['Isosurface Neighbors Grid Algorithm']);
            title("Stepped Leader Growth: Iso-grid DBM");
            App.growthAlg = 'IGDBM';
        elseif(strcmp(ev.Key, 'space'))
            runDBM();
        elseif(strcmp(ev.Key, 'backquote'))
            %deform this model
            %get rid of degenerate vertice
            dupeInds = App.AE(:, 1) == App.AE(:, 2);
            nE = App.AE(~dupeInds, :);
            deformLightningGUI(App.AV, nE, [1; size(App.AV, 1)]);
        end
    end

    function runDBM()
        if (strcmp(App.growthAlg, 'MCDBM'))
            MonteCarloDBM();
        elseif (strcmp(App.growthAlg, 'IGDBM'))
            IsoGridDBM();
        end;
    end

    function MonteCarloDBM()
        dist = 1/40;                   %distance from aggregate to sample from
        numWalks = 20;
        samplingDensity = 100;
        eps = 1/400;                 %minimum radius for WoS solve.
        
        App.not_done = true;
        while(App.not_done)
            [App.AV, App.AE, newV] = MCDBMStep(App.AV, App.AE, App.SV, App.SE, dist, numWalks, samplingDensity, eps);     
            
            [d, i, cp] = point_mesh_squared_distance(newV,App.SV, App.SE);
            if (d < dist^2)
                App.AV = [App.AV; App.SV(i, :)];
                App.AE = [App.AE; size(App.AV, 1)-1 size(App.AV, 1)];
                App.not_done = false;
            end
            
            delete(App.AEplot)
            App.AEplot = plot_edges(App.AV, App.AE, 'Color', [0.9 0.8 0.1]);
            drawnow;
        end


    end

    function  IsoGridDBM()
        res = 40 + 1;
        dist = 1/res;
        eta = 3;
        samplingDensity = 1000;
        [gridV, gridF] = create_regular_grid(res);
        L = cotmatrix(gridV, gridF);
          
        App.not_done = true;
        while(App.not_done)
            bV = [ App.SV; App.AV];
            bE = [ App.SE; App.AE+size(App.SV, 1)];
            bc = [ones(size(App.SE, 1), 1);...
                zeros(size(App.AV, 1), 1)];              %zeros at aggregate, ones at source
            gridVI = snap_points(bV, gridV);
            [gridVI, IA] = unique(gridVI);
            gridBC = bc(IA);
            
           
            randP = sampleVertIsosurface(App.AV, App.AE, dist, samplingDensity);
            [~, randPI, ~] = point_mesh_squared_distance(randP, gridV, gridF);
            randPBC = barycentric_coordinates(randP, gridV(gridF(randPI, 1), :),gridV(gridF(randPI, 2), :),gridV(gridF(randPI, 3), :));
            
            gridPhi = min_quad_with_fixed(L, [], gridVI, gridBC);
            phiRandP =  gridPhi(gridF(randPI, 1)).*randPBC(:, 1) +  gridPhi(gridF(randPI, 2)).*randPBC(:, 2) +  gridPhi(gridF(randPI, 3)).*randPBC(:, 3);
            
            
            W = phiRandP;
            W = W .^ eta;
            W = W ./ sum(W);
            %Sample one of them
            I = weightedSampling(W, 1);
            
%             tau = 1./ phiRandP;
%             [minVal, minI] = min(tau);
%             minInds = tau == minVal;
%             minVals = tau(minInds);
%             si = randi([1 length(minVals)], 1, 1);
%             I = find(minInds);
%             I = I(si);
            newV = randP(I, :);
            [AV2, AE2] = growAggregateVert(App.AV, App.AE, newV);
            App.AV = AV2; App.AE = AE2;
            
            
            [d, i, cp] = point_mesh_squared_distance(newV,App.SV, App.SE);
            if (d < dist^2)
                App.AV = [App.AV; App.SV(i, :)];
                App.AE = [App.AE; size(App.AV, 1)-1 size(App.AV, 1)];
                App.not_done = false;
            end
            
            delete(App.AEplot)
            App.AEplot = plot_edges(App.AV, App.AE, 'Color', [0.9 0.8 0.1]);
            drawnow;
        end
    end
end

