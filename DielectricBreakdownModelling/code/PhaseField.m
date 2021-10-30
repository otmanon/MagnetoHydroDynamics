% This script implements the "Phase-Field MOdel of Dielectric Breakdown in Solids",.
clear; clc;
res = 50;

[V, F] = create_regular_grid(res);

yi = (ceil(0.45*res):floor(0.55*res))';
%yi = (floor(1):floor(1))';
xi = repelem(ceil(res/2), length(yi), 1);
ii = sub2ind([res res], yi, xi);        %starting index of breakdown
phi = zeros(size(V, 1), 1); %potential
s = ones(size(V, 1), 1);    %phase field variable
s(ii) = 0;


fontsize = 24;
fig = figure('Position', [0, 0, 1500, 600]);
s1 = subplot(2, 3, [1 2]);
hold on;
tPlot = tsurf(F, V, 'CData', phi, fphong);
title("Temperature Field", 'FontSize', fontsize);
colorbar();
caxis([min(phi), max(phi)]);
xlim([-1 2]);
colormap(parula(31));

s2 = subplot(2, 3, 3);
hold on;
sPlot = tsurf(F, V, 'CData', s, fphong, falpha(1, 0));
title("Phase Field", 'FontSize', fontsize);
colorbar();
caxis([0, 1]);
colormap(parula(31));

s3 = subplot(2, 3, 4);
hold on;
colorbar();
esPlot = tsurf(F, V, 'CData', phi, fphong, falpha(1, 0));
title("ElectroStatic Stress", 'FontSize', fontsize);

s4 = subplot(2, 3, 5);
hold on;
colorbar();
iePlot = tsurf(F, V, 'CData', phi, fphong, falpha(1, 0));
title("Ionization Energy", 'FontSize', fontsize);

s5 = subplot(2, 3, 6);
hold on;
colorbar();
psPlot = tsurf(F, V, 'CData', phi, fphong, falpha(1, 0));
title("Phase Smoothing", 'FontSize', fontsize);

[EX, EY] = macGradient(reshape(phi, [res, res]), 1);

E = [EX(:), EY(:)];
    

botI = find(V(:, 2) < 0.01);
topI = find(V(:, 2) > 0.99);
BI = [botI; topI];
BPhi = [100*ones(length(botI), 1); zeros(length(topI), 1)] ;   

f = 4.*s.^3 -3.*s.^4;           %interpolation function between 0 and 1
fp = 12.*s.^2 - 12.*s.^3;
eta = 0.001;  

h = 1/(res - 1);
Ls = fdLaplacianMatrix(res, ones(length(s),1), 1);
step = 0;
eps0 = 1;
m = 1;
gamma = 1;
l = h;
eps = eps0 ./(f + eta);

noise = (rand(res)).*0.2 -  0.1;
noise = noise(:);

deltaT = 0.001;
while(true)
    
    LPhi = fdLaplacianMatrix(res, eps, h);% cotmatrix(V, F);
    phi = min_quad_with_fixed(LPhi, [], BI, BPhi);  %solve for potential field;
     
    tPlot.CData = phi;
    tPlot.Vertices = [V, phi];
 
    [EX, EY] = macGradient(reshape(phi, [res, res]), 1);
    E = [EX(:), EY(:)];
    
    ESTerm = -eps0 *fp./(2*(f + eta).^2).* (dot(E, E, 2));%-
    ITerm =  - gamma*(2.*s - 6.*s.^2 + 4.*s.^3);
    PTerm = gamma*0.5*Ls*s;
    
    esPlot.CData = ESTerm;
    esPlot.Vertices = [V, ESTerm];
    iePlot.CData = ITerm;
    iePlot.Vertices = [V, ITerm];
    psPlot.CData = PTerm;
    psPlot.Vertices = [V, PTerm];
    deltaS = ESTerm + ITerm + PTerm;
    s = s + deltaT * deltaS;
    s(ii) = 0;
    
    f = 4.*s.^3 -3.*s.^4;           %interpolation function between 0 and 1
    fp = 12.*s.^2 - 12.*s.^3;
  
    alphas = 1 ./(f + eta);
    
     %Calculate E field
   % E = fdGrad(phi, res); %alternative
    E = fd_grad([res res]);
  %  [EX, EY] = macGradient(reshape(phi, [res, res]), 1);
   % E = [EX(:), EY(:)];
    
    sPlot.CData = s;
    sPlot.Vertices = [V, s];
    drawnow;
    
    if (mod(step, 100) == 0)
      %  figgif("results/PhaseField0dot0000001.gif");
    end
    step = step + 1;
end

s2 = subplot(1, 2, 2);
hold on;
colorbar();
%L = cotmatrix(V, F);

               % boundary conditions according ot x coord


tPlot = tsurf(F, V, 'CData', phi, fphong, falpha(1, 0));



hold off