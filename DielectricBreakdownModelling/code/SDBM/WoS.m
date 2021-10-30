
function [phi, gradPhi] = WoS(P, V, E, bc, eps, maxIter)
%Performs walk on spheres to solve the laplace equation with dirichlet BC
%for points P, with boundary defined by V, and E (indices into V), and
%Boundary values defined by bc.   If a walk times out (does more than 15
%bounces) then set BC=1;

%Inputs:
%   P:      List of 2D query points
%   V:      List of 2D vertices
%   E:      List of Edges indexing V
%   bc:     list of dirichlet boundary conditions aligned with V
%   eps:    stopping threshold radius for WoS algo.
   

    lengths = edge_lengths(V, E);
    maxSamples = 60;                            %If a sphere walks more than 15 times, assume it's done.
    phiArr = zeros(size(P, 1), maxIter);
    xDirArr = zeros(size(P, 1), maxIter);
    yDirArr = zeros(size(P, 1), maxIter);

    timedOut = false(size(P, 1), maxIter);
    parfor iter = 1:maxIter
        sampledDirs = zeros(size(P));
        sampleP = P;
        N = ones(size(P, 1), 1);
        R = zeros(size(P, 1), 1);
        I = zeros(size(P, 1), 1);
        C = zeros(size(P, 1), 2);
        S = zeros(size(P, 1), 1);
        
        timedOutThisIter = false(size(P, 1), 1);
        done = false(size(P, 1), 1);
        firstJump = true;
        while(1)
            [sR, sI, sC] = point_mesh_squared_distance(sampleP(~done, :), V, E);  
            sR = sqrt(sR);
            %I is the index of the edge?

            R(~done) = sR; 
            I(~done) = sI;
            C(~done, :) = sC;

            %get value
            sS = V(E(sI, 1), :) - sC;
            sS = vecnorm(sS, 2, 2)./lengths(sI);
            S(~done) = sS;

           
            done(~done) =  sR < eps;

            %
            %Are all points sufficiently close
            if all(done)
                break;
            end

            %Sample anew
            randTheta = rand(sum(~done), 1)*2*pi;
            circleOffsets = [R(~done).*cos(randTheta),  R(~done).*sin(randTheta)];
            if (firstJump)
                sampledDirs(~done, :) = circleOffsets;
                sampledDirs(~done, :) = sampledDirs(~done, :) ./ vecnorm(sampledDirs(~done, :), 2, 2);
                firstJump = false;
            end;
            sampleP(~done, :) = sampleP(~done, :) + circleOffsets;
            N(~done) = N(~done) + 1;
            done(N > maxSamples) = 1;
            timedOutThisIter(N>maxSamples) = 1;
        end
        xDirArr(:, iter) = sampledDirs(:, 1);
        yDirArr(:, iter) = sampledDirs(:, 2);
        
        S(isnan(S)) = 0;
        phiArr(:, iter) =  (1 - S) .* bc(E(I, 1)) +  S .* bc(E(I, 2)); 
        timedOut(:, iter) = timedOutThisIter;
    end;
    
    phiArr(timedOut) = 0;  %unbounded number of bounces..
    phi = mean(phiArr, 2);
    gradXSamples = xDirArr.*phiArr;
    gradYSamples = yDirArr.*phiArr;
    dirX = mean(gradXSamples, 2);
    dirY = mean(gradYSamples, 2);
    
    gradPhi = [dirX dirY];
    
end