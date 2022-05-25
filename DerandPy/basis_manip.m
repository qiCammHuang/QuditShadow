function basis_list = basis_manip(basis_list_py)
%BASIS_MANIP Manipulate basis list
%   
global Nq

basis_list = zeros(length(basis_list_py), Nq);

for i = 1:length(basis_list_py)
    one_basis_py = basis_list_py(i);
    basis_list(i, :) = cellfun(@int64, cell(one_basis_py.pop()));
end

end
