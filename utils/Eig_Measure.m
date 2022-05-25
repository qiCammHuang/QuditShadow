function mu = Eig_Measure(basis_list, dim)
%EIG_MEASURE  outputs a realization of basis measurement
%

global Nq Deigens prob_expm_accum

mu = zeros(size(basis_list,1), Nq);
for i = 1:size(basis_list, 1)
    basis = basis_list(i,:);
    
    % the line to check in measPrescribe
    basisLine = 1;
    for jrev = 0:Nq-1
        basisLine = basisLine + (basis(Nq-jrev)-1) * (dim^2 - 1)^jrev;
    end
    outcomeList = prob_expm_accum(basisLine, :);
    
    % generate random outcome for measurement
    seed = rand(1);
    sum = 0;
    for j = 1 : dim^Nq
        sum = sum + outcomeList(j);
        if seed <= sum     
           break;
        end
    end
    outcomeStr = dec2base(j-1, dim, Nq);

    % register results in mu
    for k = 1:Nq
        mu(i,k) = Deigens(base2dec(outcomeStr(k), dim) + 1, basis(k)+1);
    end
end
end