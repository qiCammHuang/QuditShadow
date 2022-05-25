function fweight = get_fweight_random(basis_list, pauli_vec)
%GET_FWEIGHT_RANDOM   Get f-weights for LBCS
%

global Nq beta_dis

for i = 1:size(basis_list,1)
    %%%%%%%%%%% for j = 1:num_observ
    for j = 1:size(pauli_vec,1)
        fi = 1;
        for loc = 1:Nq
            if pauli_vec(j,loc) == 0
                % fi = fi*1;
            elseif pauli_vec(j,loc) == basis_list(i,loc)
                fi = fi/beta_dis(loc,basis_list(i,loc));
            else
                fi = fi*0;
                break;
            end  
        end
        fweight(i,j) = fi; % j - > Q
    end
end
end

