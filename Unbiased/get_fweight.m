function fweight = get_fweight(basis_list,pauli_vec, dim)
%% GET_FWEIGHT  Get f-weights for Unbiased Shadow
%

global Nq
for i = 1:size(basis_list,1)
    %%%%%%%%%%% for j = 1:num_observ
    for j = 1:size(pauli_vec,1)
        fi = 1;
        for loc = 1:Nq
            if pauli_vec(j,loc) == 0
                fi = fi;
            elseif pauli_vec(j,loc) == basis_list(i,loc)     
                if nargin < 3
                    fi = fi*3;
                else
                    fi = fi*(dim^2-1);
                end
            else
                fi = fi*0;
                break;
            end
            
        end
        fweight(i,j) = fi;
    end
end
end