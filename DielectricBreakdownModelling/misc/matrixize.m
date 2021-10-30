function mat=matrixize(vec, res)
%Given a flattened vector (previously flattened by vectorize), unflatten it
    mat = zeros(res, res);
    mat(:) = vec;
    mat = mat';    
end