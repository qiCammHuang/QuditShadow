function psi = get_GHZ()

global Nq dim

% Initialise psi
psi = zeros(dim^Nq, 1);

% |psi> = |GHZ>
for i = 1:dim
    string = '';
    for j = 1:Nq
        string = [string, num2str(i-1)];
    end
    index = base2dec(string, dim);
    psi(index + 1) = 1;
end

% Normalise
psi = sparse(psi)/norm(psi);
end