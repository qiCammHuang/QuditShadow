clear all

alllines = 30;

hamil = int64(zeros(alllines, 5));

cnt_line = 1;

while cnt_line < alllines + 1
    towrite = int64(zeros(1,5));

    a = randi([1,4]);
    b = randi([1,4]);
    
    if a == b
        h = randi([1,3]);
        towrite(1, a) = h;
        towrite(1, 5) = int64(1);
    else
        h1 = randi([1,3]);
        h2 = randi([1,3]);

        towrite(1, a) = h1;
        towrite(1, b) = h2;
        towrite(1, 5) = int64(1);
    end

    flag = 0;
    for i = 1:size(hamil,1)
        if towrite == hamil(i,:)
            flag = 1;
            break
        end
    end

    if flag == 0
        hamil(cnt_line, :) = towrite;
        cnt_line = cnt_line + 1;
    end
end

writematrix(hamil, 'dHamil.txt', 'Delimiter', 'tab');
type 'dHamil.txt';



