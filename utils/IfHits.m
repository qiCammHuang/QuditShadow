function res = IfHits(observ, basis)
%IFHITS Shows if basis hits observable
%   
global Nq

res = 1;
for i = 1:Nq
    if (observ(i) > 0.5) && (abs(observ(i) - basis(i)) > 0.5)
        % observ(i) is not identity, while basis(i) is not observ(i)
        % does NOT hit!
        res = 0;
        break
    end
end

end

