function beta_py = bias_manip(beta_matlab)
%BIAS_MANIP Manipulate beta distribution
%
%   Note that we intentionally add a bias of 0.0 to identity

beta_py = py.list();

for i = 1:size(beta_matlab, 1)
    one_beta_py = py.list(double([0.0, beta_matlab(i,:)]));
    beta_py.append(one_beta_py);
end

end

