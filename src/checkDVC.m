function input = checkDVC(input, x, isProposed)

[tmp strNum] = size(input);

if ~exist('A.lt') % for Merrit
    A.lt = [];
end
if ~exist('A.ut') % for Merrit
    A.ut = [];
end
fprintf('== DVC Information ==\n');

for str = 1:strNum
    A = input{str};
    dose = A.mat * x;
    numVox = A.size;
    input{str}.le = zeros(size(A.lp));
    input{str}.ue = zeros(size(A.up));
    doseSort = sort(dose,'descend');
    for ind = 1 : max(size(A.lp))
        topNumVox = floor(input{str}.lp(ind) * numVox);
        input{str}.lDVC(ind) = doseSort(topNumVox);
        fprintf('L_%d^%.2f = %.1f  => %.1f : ', str, input{str}.lp(ind), ...
            input{str}.ld(ind), input{str}.lDVC(ind));
        if input{str}.lDVC(ind) >= input{str}.ld(ind) - (1.0e-2)
            fprintf(' Pass ');
        else
            fprintf(' Fail ');
        end
        fprintf('\n');
    end

    
    
    for ind = 1 : max(size(A.up))
        topNumVox = ceil(input{str}.up(ind) * numVox);
        input{str}.uDVC(ind) = doseSort(topNumVox);
        fprintf('U_%d^%.2f = %.1f  => %.1f : ', str, input{str}.up(ind), ...
            input{str}.ud(ind), input{str}.uDVC(ind));
        if input{str}.uDVC(ind) <= input{str}.ud(ind) + (1.0e-2)
            fprintf(' Pass ');
        else
            fprintf(' Fail ');
        end
        fprintf('\n');
    end
end

if isProposed ~= 1
    return;
end

fprintf('== Hot Spot Information ==\n');
for str = 1:strNum
    A = input{str};
    dose = A.mat * x;
    numVox = A.size;
    for ind = 1 : max(size(A.lp))
        input{str}.lHot(ind) = length(find(dose < input{str}.lt(ind)));
        fprintf('L_%d^%.2f = %.1f  => %.1f : ', str, input{str}.lp(ind), ...
            input{str}.ld(ind), input{str}.lt(ind));
        fprintf('Num = %d/%d, Per = %.2f%%\n', input{str}.lHot(ind), ...
            numVox, input{str}.lHot(ind) / numVox*100);
    end
    for ind = 1 : max(size(A.up))
        input{str}.uHot(ind) = length(find(dose > input{str}.ut(ind)));
        fprintf('U_%d^%.2f = %.1f  => %.1f : ', str, input{str}.up(ind), ...
            input{str}.ud(ind), input{str}.ut(ind));
        fprintf('Num = %d/%d, Per = %.2f%%\n', input{str}.uHot(ind), ...
            numVox, input{str}.uHot(ind) / numVox*100);
    end
end



