%Code taken from Seungbae Bang, git repo - vector_graphics_on_sphere
function [BV,BE,BN] = linspaced_square(O,r,samples_per_edge)

    %O = [0 0];
    %r = 2;

    C = [O(:,1)+r,O(:,2)-r; O(:,1)+r,O(:,2)+r; ...
         O(:,1)-r,O(:,2)+r; O(:,1)-r,O(:,2)-r];

    E = [1,2;2,3;3,4;4,1];

    %samples_per_edge = 20;
    t = linspace(0.0,1.0,samples_per_edge+1);
    t = t(1:(end-1));
    t = reshape(repmat(t,2,1),1,size(t,2)*2);
    t = repmat(t,size(E,1),1);

    sp = repmat(C(E(:,1),:),1,samples_per_edge);
    ep = repmat(C(E(:,2),:),1,samples_per_edge);
    S = sp.*(1-t) + ep.*t;
    BV = reshape(S',2,prod(size(S'))/2)';
    id = 1:size(BV,1);
    nid = circshift(id,-1);
    BE = [id; nid]';

    A = BV(BE(:,2),:)-BV(BE(:,1),:);
    UN = A*[0 -1;1 0];
    BN = normalizerow(UN);

end