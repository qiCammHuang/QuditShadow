function [En, Ol] = postShadow_Hamil(mu,Sample,list,fweight)
%POSTSHADOW_HAMIL   Generate estimate of energy and each operator
%
%   mu: raw measurement outcomes
%   Sample: number of samples
%   list: a list, 1:num_observ
%   fweight: corresponding shadow f-weights
%

global pauli_vec

Estimate = zeros(Sample, length(list));

parfor i = 1:Sample
    Energy(i) = 0;

    for j = list
        % for the j-th observable
        mu_shadow = 1;
        vec = pauli_vec(j,1:end-1);

        index = find(vec~=0);
        for tmp = 1:length(index)
            mu_shadow = mu_shadow * mu(i,index(tmp));
        end
        
        Estimate(i,j) = fweight(i,j) * mu_shadow;
        Energy(i) = Energy(i) + pauli_vec(j,end)*Estimate(i,j);
    end
end

En = 0;
for i = 1:Sample
    En = En + real(Energy(i));
end
En = En/Sample;

if nargout == 2
    Ol = real(mean(Estimate, 1));
end

end
