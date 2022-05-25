function [observ_py, weight_py] = hamil_manip(observ_matlab)
%HAMIL_MANIP    Manipulate Hamiltonian format
%

global Nq

observ_py = py.list();
weight_matlab = [];

for i = 1:size(observ_matlab,1)
    one_observ_py = py.list();
    for j = 1:Nq
        if observ_matlab(i,j) > 0.5
            one_observ_matlab = [int64(observ_matlab(i,j)), int64(j-1)];
            one_observ_py.append(py.tuple(py.list(one_observ_matlab)));
        end
    end
    observ_py.append(one_observ_py);
    weight_matlab = [weight_matlab, abs(observ_matlab(Nq+1))];
end
weight_matlab = weight_matlab / max(weight_matlab);

weight_py = py.list(weight_matlab);
end