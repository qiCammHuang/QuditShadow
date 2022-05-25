function basis_list = get_basis_list(Sample, dim)
%BASIS_LIST   Get LBCS List of Bases
%   

global Nq beta_dis

basis_list = zeros(Sample,Nq);
for i = 1:Sample
    for j = 1:Nq
        q = rand(1);
        if nargin == 1
            if q < beta_dis(j,1)
                basis_list(i,j) = 1;
            elseif q < beta_dis(j,1)+beta_dis(j,2)
                basis_list(i,j) = 2;
            elseif q < beta_dis(j,1)+beta_dis(j,2)+beta_dis(j,3)
                basis_list(i,j) = 3;
            else
                'Error in basis list!!'
            end
        else
            sum = 0;
            for k = 1:dim^2-1
                sum = sum + beta_dis(j,k);
                if q <= sum
                    basis_list(i,j) = k;
                    break                 
                end
            end
        end
             
    end
end
end