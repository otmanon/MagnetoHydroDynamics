function vec=vectorize(mat)
%flattens matrix to vector in row major order
    tempMat = mat';
    vec = (tempMat(:));
    
end