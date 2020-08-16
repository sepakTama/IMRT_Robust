function [Input,x] = lambdaIteration(Input, alpha, iterationNum, lambdaIte, option, priority)

% lambdaIte ... the number of iteration in order to deciding lambda
% option: To pass the option to the next function(PenaltyIterate.m)
% priority: if priority is 'non', then all PTV and OAR are equal to handle
%          if priority is 0, then don't update the parameter


[~, strNum] = size(Input);

% conv: converge rate in each iteration
% stepsize: the size of prameter in each iteration
% rate: the ratio to step size powed by "priority (Num)"
% endoption: if endoption is 0, it breaks the iteration when all deviations are nonpositive
eps = 1e-3;
conv = 0.95;
rate = 0.1;
stepsize = 100;
% 臓器ごとにstepsizeを変えるならば、stepsizeを臓器分配列の形で用意。
endOption = 1;


if priority == 0
    rate = 1;
end

lambdaMinus = [];
lambdaPlus = [];

% set parameters directory
for str = 1:strNum
    Input{str}.LMinus = [];
    Input{str}.LPlus = [];
    fprintf('===Structure{%d}===\n', str);   
    for ind = 1:max(size(Input{str}.lp))
        fprintf('->Lp = {%.2f}\n', Input{str}.lp(ind));
        prompt = 'Enter the parameters of lambdaMinus.\n';
        Input{str}.LMinus = [Input{str}.LMinus; input(prompt)];
    end
    for ind = 1:max(size(Input{str}.up))
        fprintf('->Up = {%.2f}\n', Input{str}.up(ind));
        prompt = 'Enter the parameters of lambdaPlus.\n';
        Input{str}.LPlus = [Input{str}.LPlus; input(prompt)];
    end
    % print
    lambdaMinus = [lambdaMinus; Input{str}.LMinus];
    lambdaPlus = [lambdaPlus; Input{str}.LPlus];
end

% log the data of paramter
fileID1 = fopen('./track/trash.txt','a');
fileID2 = fopen('./track/trash_data.txt','a');

data_Plus_x = [];
data_Minus_x = [];
data_Plus_y = [];
data_Minus_y = [];


lambda_total = tic;


for ite = 1:lambdaIte
    % solve the optimization problem
    [Input, x, opt,  t_history, ~, ~] = PenaltyIterate(Input, alpha, iterationNum, option, lambdaMinus, lambdaPlus);
    
    fprintf('===lambdaIteration{%d}===\n', ite);
    % get the gap between DVC and optimal value
    [lowGap, highGap] = Compare(Input, opt);
    
    lambda_time = toc(lambda_total);
   
    % log parameters, time about n times optimization prob., gap between DVC and optim value, and lambda_time
    fprintf(fileID1, '=======lambdaite %d=====\r\n', ite);
    fprintf(fileID1, '/lambdaPlus\n');fprintf(fileID1,'%f\n', lambdaPlus);
    fprintf(fileID1, '/lambdaMinus\n');fprintf(fileID1,'%f\n', lambdaMinus);
    fprintf(fileID1, '/t_hisotry\n');fprintf(fileID1,'%f\n', t_history);
    fprintf(fileID1, '/lowGap\n');fprintf(fileID1,'%f\n', lowGap);
    fprintf(fileID1, '/highGap\n');fprintf(fileID1,'%f\n', highGap);
    fprintf(fileID1, '/lambda_time\n');fprintf(fileID1,'%f\n',lambda_time);

    data_Plus_x = [data_Plus_x; lambdaPlus];
    data_Minus_x = [data_Minus_x; lambdaMinus];
    data_Plus_y = [data_Plus_y; highGap];
    data_Minus_y = [data_Minus_y; lowGap];
    

    % if satisfy all gap are negative, then stop iteration
    isOK = endOption;
    for i = 1:size(lowGap) % count sum of the number of deviation 
        if lowGap(i) > 0
            isOK = isOK + 1;
            break;
        end
    end
    for i = 1:size(highGap)
        if highGap(i) > 0
            isOK = isOK + 1;
            break;
        end
    end
    
    if isOK == 0 
        break;
    else
        lambdaMinus = [];
        lambdaPlus = [];
        UIndex = 0;
        LIndex = 0;
        if priority == 'non'
            for str = 1:strNum
                for ind = 1:max(size(Input{str}.lp))
                    if lowGap(LIndex + ind) < 0
                        Input{str}.LMinus(ind) = max(0,Input{str}.LMinus(ind) - stepsize);
                    elseif lowGap(LIndex + ind) > 0
                        Input{str}.LMinus(ind) = Input{str}.LMinus(ind) + stepsize;
                    end
                    lambdaMinus = [lambdaMinus; Input{str}.LMinus(ind)];
                end

                for ind = 1:max(size(Input{str}.up))
                    if highGap(UIndex + ind) < 0
                        Input{str}.LPlus(ind) = max(0,Input{str}.LPlus(ind) - stepsize);
                    elseif highGap(UIndex + ind) > 0
                        Input{str}.LPlus(ind) = Input{str}.LPlus(ind) + stepsize;
                    end
                    lambdaPlus = [lambdaPlus; Input{str}.LPlus(ind)];
                end
            LIndex = LIndex + max(size(Input{str}.lp));
            UIndex = UIndex + max(size(Input{str}.up));
            end    
        else
            for str = 1:strNum
                if priority(str) > 0
                    for ind = 1:max(size(Input{str}.lp))
                        if lowGap(LIndex + ind) < 0
                            Input{str}.LMinus(ind) = max(0,Input{str}.LMinus(ind) - (stepsize*rate.^(priority(str)-1))); %割合をpriorityの乗数にした。
                        elseif lowGap(LIndex + ind) > 0
                            Input{str}.LMinus(ind) = Input{str}.LMinus(ind) + (stepsize*rate.^(priority(str)-1));
                        end
                        lambdaMinus = [lambdaMinus; Input{str}.LMinus(ind)];
                    end
                elseif priority(str) < 0 % 健全な臓器
                    for ind = 1:max(size(Input{str}.lp))
                        Input{str}.LMinus(ind) = Input{str}.LMinus(ind) + (stepsize*rate.^(abs(priority(str))-1));
                        lambdaMinus = [lambdaMinus; Input{str}.LMinus(ind)];
                    end
                else
                    for ind = 1:max(size(Input{str}.lp))
                        lambdaMinus = [lambdaMinus; Input{str}.LMinus(ind)];
                    end
                end
                
                if priority(str) > 0
                    for ind = 1:max(size(Input{str}.up))
                    
                        if highGap(UIndex + ind) < 0
                            Input{str}.LPlus(ind) = max(0,Input{str}.LPlus(ind) - (stepsize*rate.^(priority(str)-1)));
                        elseif highGap(UIndex + ind) > 0
                            Input{str}.LPlus(ind) = Input{str}.LPlus(ind) + (stepsize*rate.^(priority(str)-1));
                        end
                        lambdaPlus = [lambdaPlus; Input{str}.LPlus(ind)];
                    end
                elseif priority(str) < 0 %健全な臓器
                    for ind = 1:max(size(Input{str}.up))
                        Input{str}.LPlus(ind) = Input{str}.LPlus(ind) + (stepsize*rate.^(abs(priority(str))-1));
                        lambdaPlus = [lambdaPlus; Input{str}.LPlus(ind)];
                    end
                else
                    for ind = 1:max(size(Input{str}.up))
                        lambdaPlus = [lambdaPlus; Input{str}.LPlus(ind)];
                    end
                end
            LIndex = LIndex + max(size(Input{str}.lp));
            UIndex = UIndex + max(size(Input{str}.up));
            end
        end

    end
        %stepsizeを収束に持っていくためにstepsizeを小さくする。
    stepsize = stepsize * conv;   
end
    lambda_time = toc(lambda_total);
    
    fprintf(fileID2, 'data_Plus_x = ['); fprintf(fileID2, '%f, ', data_Plus_x); fprintf(fileID2,'];'); fprintf(fileID2,'\n'); 
    fprintf(fileID2, 'data_Minus_x = ['); fprintf(fileID2, '%f, ', data_Minus_x); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    fprintf(fileID2, 'data_Plus_y = ['); fprintf(fileID2, '%f, ', data_Plus_y); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    fprintf(fileID2, 'data_Minus_y = ['); fprintf(fileID2, '%f, ', data_Minus_y); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    
    fprintf('*** %d lambda Iteration ***\n', ite);
    lambda_time
end
