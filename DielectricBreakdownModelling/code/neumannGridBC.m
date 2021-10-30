function [Ae, B] = neumannGridBC(res, h, bc)
%Gives Ae matrix for linear equality constraints that correspond to zero
%neumann

    ijv = [];
    
%     % bc neumann on top
%     tyi = repelem([res], res-2, 1 );
%     txi = (2:res-1)';
%     ti = sub2ind([res, res], tyi, txi);
%     
%     tyio = repelem([res-1], res-2, 1 );
%     txio = (2:res-1)';
%     tio = sub2ind([res, res], tyio, txio);
%     
%         
%     % bc neumann on bottom
%     byi = repelem(1, res-2, 1 );
%     bxi = (2:res-1)';
%     bi = sub2ind([res, res], byi, bxi);
%     
%     byio = repelem(2, res-2, 1 );
%     bxio = (2:res-1)';
%     bio = sub2ind([res, res], byio, bxio);
    
    
    %bc neumann left
    lxi = repelem(1, res, 1 );
    lyi = (1:res)';
    li = sub2ind([res, res], lyi, lxi);
    
    lxio = repelem(2, res, 1 );
    lyio = (1:res)';
    lio = sub2ind([res, res], lyio, lxio);
    
    
    %bc neumann right
    rxi = repelem(res, res, 1 );
    ryi = (1:res)';
    ri = sub2ind([res, res], ryi, rxi);
    
    rxio = repelem(res-1, res, 1 );
    ryio = (1:res)';
    rio = sub2ind([res, res], ryio, rxio);
    
%     bdio = [tio; bio;lio;rio];
%     bdi = [ti; bi; li; ri];

    bdio = [lio;rio];
   
    bdi = [li; ri];

    numConstraints = length(bdio);
    
    nbc = zeros(numConstraints,1);
    nbc(numConstraints/2 + 1:numConstraints) = -bc;
    nbc(1:numConstraints/2) = bc;
    
    ijv = [(1:numConstraints)', bdi,  ones(numConstraints, 1)];
    ijv = [ijv;
            (1:numConstraints)', bdio, -ones(numConstraints, 1)];
    
%     numConstraints = numConstraints + 1;
%     ijv = [ijv; numConstraints, res*res/2 + floor(res/4), 1];
%  %   numConstraints = numConstraints + 1;
 %   ijv = [ijv; numConstraints, res*res/2 + floor(res*0.75), 1];
    Ae = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3), numConstraints, res*res);
    B = zeros(numConstraints, 1);
    B = nbc;
  
    
end

