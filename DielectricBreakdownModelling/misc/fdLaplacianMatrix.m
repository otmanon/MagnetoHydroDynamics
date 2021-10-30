function L = fdLaplacianMatrix(res, alphas, h)
%FDLAPLACIAN Summary of this function goes here
%   Calculates the finite difference laplacian matrix for a unit grid with
%   res vertices per row/col. Assumes 0-Neumann BC

    L = sparse((res+2)*(res + 2));
        
    if (length(alphas) < (res+2)*(res + 2))
        
        sind = 1:length(alphas);
        [suby, subx] = ind2sub(res, sind);
        ralphas = zeros((res+2) * (res+2));
        rind = sub2ind([res+2 res+2], suby+1, subx+1);
        ralphas(rind) = alphas;
        alphas = ralphas;
    end;
    
    ijv = [];
    indeces=[];
    aci = 1;
    for x=2:res + 1
        for y=2:res +1
            index = sub2ind([res+2 res+2], y, x);
            indeces = [indeces; index];
           
            
            ti = sub2ind([res+2 res+2], y+1, x);
            at = 0.5*(alphas(index) + alphas(ti));
            
            bi = sub2ind([res+2 res+2], y-1, x);
            ab = 0.5*(alphas(index) + alphas(bi));
            
            ri = sub2ind([res+2 res+2], y, x+1);
            ar = 0.5*(alphas(index) + alphas(ri));
            
            li = sub2ind([res+2 res+2], y, x-1);
            al = 0.5*(alphas(index) + alphas(li));
      
            ijv = [ijv;...
                ti ti -at;ti index at;...
                bi bi -ab;bi index ab;...
                ri ri -ar;ri index ar;...
                li li -al;li index al];
                
            
   
                   
        end
    end
    L = sparse(ijv(:, 1), ijv(:, 2), ijv(:, 3));
    L = L(indeces, indeces)/(h*h);
  

end



