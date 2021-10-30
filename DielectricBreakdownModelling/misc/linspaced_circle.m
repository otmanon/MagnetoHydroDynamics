%Code taken from Seungbae Bang, git repo - vector_graphics_on_sphere
function [C,E,N] = linspaced_circle(O,r,ns)

    theta = linspace(0,2*pi,ns+1);theta = theta(1:end-1)';
    C=[cos(theta) sin(theta)];
    C=r*C;
    C=C+repmat(O,size(C,1),1);

    id = 1:size(C,1);
    E = [id; circshift(id,-1)]';

    A = C(E(:,2),:)-C(E(:,1),:);
    UN = A*[0 -1;1 0];
    N = normalizerow(UN);

end