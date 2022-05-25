function v = get_product_state(basis_index, outcome_index)
%GET_PRODUCT_STATE  specific product state vector
%   product state of specific projection outcome, prescribed by specific
%   GGM bases and measurement outcomes
%
%   v is a column vector

global Nq dim Veigens

basisStr = dec2base(basis_index, dim^2-1, Nq);
outcomeStr = dec2base(outcome_index, dim, Nq);

v = eye(1);
for n = 1:Nq
    % Veigens(:, base2dec(outcomeStr(n), dim) + 1, base2dec(basisStr(n),dim^2-1) + 2)
    v = sparse(kron(v, sparse(Veigens(:, base2dec(outcomeStr(n), dim) + 1, base2dec(basisStr(n),dim^2-1) + 2))));
end

