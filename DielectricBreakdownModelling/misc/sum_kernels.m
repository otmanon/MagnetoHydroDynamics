function T = sum_kernels(qV, AV, SV);
%SUM_KERNELS Summary of this function goes here
%   Detailed explanation goes here
res = sqrt(size(qV, 1));
h = 1/res;
r = 1/h;

DA = pdist2(qV, AV);
ratio = DA ./ r;
ratio(DA == 0) = 0;         % be careful, if the query and aggregate point are the same, set that to zero

AT = -(1 - ratio);
AT = sum(AT, 2);

ST = zeros(size(AT, 1), 1);
if (size(SV, 1)>0)
    DS = pdist2(qV, SV);
    ratio = DS ./ r;
    ratio(DS == 0) = 0;
    ST = 1 - ratio;
    ST = sum(ST, 2);
end
T = AT + ST;

T = T - min(T);
T = T./max(T);


    
end

