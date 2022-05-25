clear all

hamil = load('cHamil.txt');

flag = 0;

for i = 1:size(hamil,1)
    for j = i+1:size(hamil,1)
        if hamil(i,:) == hamil(j,:)
            flag = 1
            i
            j
            break
            break
        end
    end
end
