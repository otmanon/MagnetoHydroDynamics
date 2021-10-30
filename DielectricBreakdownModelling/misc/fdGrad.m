function G = fdGrad(U,res)
%FDGRAD Summary of this function goes here. Uses finite differences to get
%gradient of U vector, where U is a flattened grid vector
%   Detailed explanation goes here
    Umat = reshape(U, [res, res]);
    [GX, GY] = gradient(Umat, 1/(res-1));
    G = [GX(:), GY(:)];
    
end

