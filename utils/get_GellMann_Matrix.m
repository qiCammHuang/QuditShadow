function [GGM, v, d] = get_GellMann_Matrix(dim)
%GET_GELLMANN_MATRIX    Generalised Gell-Mann Matrices
%   Return Generalised Gell-Mann Matrices, along with respecive
%   eigenvectors and eigenvalues.
%   
%   GGM(:,:,g) is the g-th generalised Gell-Mann Matrix
%   
%   with descending order wrt eigenvalues:
%   v(:,i,g) is the i-th eigenvector of the g-th GGM
%   d(i,g) is the i-th eigenvalue of the g-th GGM
%

basis = eye(dim);
GGM(:,:,1) = basis;
v = zeros(dim,dim,dim^2);
d = zeros(dim,dim^2);

offset = 2;
% off-diagonal symmetric
for k = 1:dim
    for j = 1:k-1
        GGM(:,:,offset) = basis(:,j)*basis(:,k)' + basis(:,k)*basis(:,j)';
        offset = offset + 1;
    end
end

% off-diagonal anti-symmetric
for k = 1:dim
    for j = 1:k-1
        GGM(:,:,offset) = -1j*basis(:,j)*basis(:,k)' + 1j*basis(:,k)*basis(:,j)';
        offset = offset + 1;
    end
end

% diagonal
for l = 1:dim-1
    tmp = 0;
    for j = 1:l
         tmp = tmp +  basis(:,j)*basis(:,j)';
    end
    GGM(:,:,offset) = sqrt(2/(l*(l+1))) * (tmp-l*basis(:,l+1)*basis(:,l+1)');
    offset = offset + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigen
for g = 1:dim^2
    ggm = GGM(:,:,g);

    [v_tmp, d_tmp] = eig(ggm);
    [~, ind] = sort(diag(d_tmp), "descend");

    v(:,:,g) = v_tmp(:, ind);
    d(:,g) = diag(d_tmp(ind,ind));
end

