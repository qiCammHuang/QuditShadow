function [En, Ol] = postDerand_Hamil(mu,Sample,list,basis_list)
%POSTDERAND_HAMIL   Generate estimate of energy and each operator
%
%   mu: raw measurement outcomes
%   Sample: number of samples
%   list: a list, 1:num_observ
%   basis_list: list of bases
%

global pauli_vec

Estimate = zeros(Sample, length(list));
Validity = zeros(Sample, length(list));

for i = 1:Sample
    for j = list
        % for the j-th observable
        vec = pauli_vec(j,1:end-1)
        if IfHits(vec, basis_list(i,:)) > 0.5
            % if sample hits observable, it provides actionable info
            Validity(i,j) = 1;
            mu_sift = 1;

            index = find(vec~=0);
            for tmp = 1:length(index)
                mu_sift = mu_sift * mu(i,index(tmp));
            end

            Estimate(i,j) = mu_sift;
        end % if sample doesn't hit, do nothing to this observable this round.
    end
end

if nargout == 2
    overall_valid = sum(Validity, 1);
    for l = 1:length(overall_valid)
        if overall_valid(l) < 0.5
            print('The', l, 'th observable was Not hit!')
            overall_valid(l) = overall_valid(l) + 1;
        end
    end
    Ol = real(sum(Estimate, 1) ./ overall_valid);
end

En = 0;
for l = list
    En = En + pauli_vec(l, end) * Ol(l);
end


end
