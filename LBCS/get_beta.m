function beta_dis = get_beta(Iter_num, delta, dim)

%%%% Minimise variance to get Beta distribution %%%%

global Nq pauli_vec

if nargin == 2
    beta_tmp = 1/3*ones(Nq,3);
    dim = 2;
else
    beta_tmp = 1/(dim^2-1)*ones(Nq,dim^2-1);
end

% delta = 0.1;

%%%% use convex function's optimisation %%%%
for iter = 1:Iter_num
    beta_CP = get_beta_CP(beta_tmp,dim);
    beta_tmp = (1-delta)*beta_tmp + delta*beta_CP;
    
    for i = 1:size(beta_tmp,1)
        if abs(sum(beta_tmp(i,:))-1)>1e-5
            iter 
            beta_tmp(i,:) = beta_tmp(i,:)/sum(beta_tmp(i,:));
        end
    end
    if max(beta_tmp - beta_CP) < 0.0001
        iter
        break
    end
    
end

beta_dis = beta_tmp;

end

function beta_CP = get_beta_CP(beta_dis,dim)

global Nq pauli_vec
%%%%%%%%%% beta_CP %%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nq % for every site     
    for index = 1:dim^2-1 % for every Pi
        b1 = 0;
        b2 = 0;
        % sum over Q
        for num = 1:size(pauli_vec,1) % number of observable
            if pauli_vec(num,i) == index % Q_i = P_i         
                %%%%%%%%%% Supp(Q) %%%%%%
                vec = pauli_vec(num,1:end-1);
                vec(find(vec==0))=[];
                prod = 1;
                for tmp = 1:size(vec,1)
                    prod = prod /beta_dis(tmp,vec(tmp));

                end
                b1 = b1 + pauli_vec(num,end)^2 * prod;
            end
            
            if pauli_vec(num,i) ~= 0 % Q_i ?? I        
                %%%%%%%%%% Supp(Q) %%%%%%
                vec = pauli_vec(num,1:end-1);
                vec(find(vec==0))=[];
                prod = 1;
                for tmp = 1:size(vec,1)
                    prod = prod * 1/beta_dis(tmp,vec(tmp)); 

                end
                b2 = b2 + pauli_vec(num,end)^2 * prod;
            end
  
            
        end
        beta_CP(i,index) = b1/b2; % beta_i(Pi)
    end
end

end