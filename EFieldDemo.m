function EFieldDemo()
%EFIELDDEMO Summary of this function goes here
%   Detailed explanation goes here
close all; clear all; clc;
numP = 1;

SimOBJ.PX =  rand(numP, 3);             %ChargePositions
SimOBJ.PV = rand(numP, 3) - [0.5, 0.5, 0.5];
SimOBJ.PV(:, 3) = 0;
SimOBJ.Q = randi(3, numP, 1) - 2;          %negative 1 if electron, 1 if pos ion, 0 if neutral                       %Charges
SimOBJ.Q = 1;
SimOBJ.M = ones(numP, 1);
SimOBJ.invM = 1./SimOBJ.M;
SimOBJ.dt = 5;
SimOBJ.E = [0, 0, 0];           %E field
SimOBJ.B = [0, 0, 1];           % B field

hold on;
SimOBJ.particleVis = scatter(SimOBJ.PX(:, 1), SimOBJ.PX(:, 2), 'filled', 'CData', SimOBJ.Q);
axis equal;
%     xlim([0, 1]); ylim([0, 1]);
colorbar();

title("Magnetic Field Pointing Out of the page");


mat = eye(3) - SimOBJ.dt .* SimOBJ.Q .* SimOBJ.invM .* skew(SimOBJ.B);

matZero = zeros(6);
matZero(1:3, 1:3) = eye(3)
matZero(4:6, 4:6) = mat
matZero(1:3, 4:6) = SimOBJ.dt.*mat
det(matZero)
%  set(gca,'ButtonDownFcn',@mouse_click_callback);
velVis = [];
waitforbuttonpress();
while(true)
    
    lorentzF =  (SimOBJ.E   + cross(SimOBJ.PV, SimOBJ.B)).*SimOBJ.Q;
    SimOBJ.PV = SimOBJ.dt .* lorentzF .* SimOBJ.invM + SimOBJ.PV;
    SimOBJ.PX = SimOBJ.PV .* SimOBJ.dt + SimOBJ.PX;
%     v2 = SimOBJ.PV/mat;
%     SimOBJ.PV = v2;
    
    
    scatter(SimOBJ.PX(:, 1), SimOBJ.PX(:, 2), 'MarkerFaceColor', [0.5, 0.5, 0.99], 'MarkerEdgeColor', 'None');
    
    delete(SimOBJ.particleVis);
    SimOBJ.particleVis = scatter(SimOBJ.PX(:, 1), SimOBJ.PX(:, 2), 'filled', 'CData', [1 0.2 0.2]);
    quiver(SimOBJ.PX(:, 1), SimOBJ.PX(:, 2), SimOBJ.dt*SimOBJ.PV(:, 1), SimOBJ.dt*SimOBJ.PV(:, 2), 'Color', [0 0 0]);
    %        SimOBJ.particleVis.XData = SimOBJ.PX(:, 1);
    %        SimOBJ.particleVis.YData = SimOBJ.PX(:, 2);
%     figgif("MagneticImplicitDemo.gif");
    pause(0.1);
end


    function X = skew(x)
        X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
    end
end

