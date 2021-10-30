function W = compute_skinning_weight(V,F,b)

% get constraints
bc = eye(length(b));

% get bilaplacian
L = cotmatrix(V,F);
M = massmatrix(V,F);
invM = diag(diag(M).^(-1));
Q = L * invM * L;

% minimize quadratic problem
% min 0.5*x'*Q*x - x'*B
unknown = [1:size(V,1)]';
unknown(b) = [];
LHS = Q(unknown, unknown);
RHS = -Q(unknown, b) * bc;
W_unknown = LHS \ RHS;

% assemble output 
W = zeros(size(V,1), length(b));
W(unknown,:) = W_unknown;
W(b,:) = bc;

% make sure sum up to one
W = W./repmat(sum(W,2),1,size(W,2));
