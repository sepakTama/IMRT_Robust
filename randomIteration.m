function[Input, x, optPlus, optMinus, lambda_time] = randomIteration(Input, alpha, iterationNum, lambdaIte, option, priority)

%lambdaIte ... the number of iteration to decide lambda
%rate:priorityでべき乗になる
%stepsize:各反復でのパラメータの更新サイズ
%conv:買う反復でstepsizeの収束
%endoption:0なら、すべてが負なら終了する条件を付けられる
%priority:'non'なら優先順位はない、0なら更新しない

% 設定
conv = 0.96;
rate = 0.1;
stepsize = 100;

[~, strNum] = size(Input);

eps = 1e-3;
% 臓器ごとにstepsizeを変えるならば、stepsizeを臓器分配列の形で用意。
endOption = 1;
% 0なら終了条件あり

if priority == 0
    rate = 1;
end

global init;
global init_lambdaMinus;
global init_lambdaPlus;

% set the parameter of all DVCs
if isempty(init) == 1
    lambdaMinus = [];
    lambdaPlus = [];
    
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
        lambdaMinus = [lambdaMinus; Input{str}.LMinus];
        lambdaPlus = [lambdaPlus; Input{str}.LPlus];
    end
    init_lambdaMinus = lambdaMinus;
    init_lambdaPlus = lambdaPlus;
    init = 0;    
    for i = 1:strNum % 別の\alphaを試すときはデータを再設定する必要あり    
        Input{i}.or_mat = Input{i}.mat;
        Input{i}.or_val = Input{i}.val;
    end
elseif isempty(init) == 0
    % 初期化
    lambdaMinus = init_lambdaMinus;
    lambdaPlus = init_lambdaPlus;
    for str = 1:strNum
        Input{str}.LMinus = init_lambdaMinus(str);
        Input{str}.LPlus = init_lambdaPlus(str);
    end
end

fileID1 = fopen('trash.txt','a');
fileID2 = fopen('trash_data.txt','a');


data_Plus_x = [];
data_Minus_x = [];
data_Plus_y = [];
data_Minus_y = [];

lambda_total = tic;

for ite = 1:lambdaIte
    for i = 1:strNum % データの修復
        Input{i}.mat = Input{i}.or_mat;
        Input{i}.val = Input{i}.or_val;
    end        
    [Input, x, opt,  t_history, ~, ~] = PenaltyIterate(Input, alpha, iterationNum, option, lambdaMinus, lambdaPlus);
    fprintf('===lambdaIteration{%d}===\n', ite);
    [lowGap, highGap] = Compare(Input, opt);
    
    lambda_time = toc(lambda_total);
   
    % write the data
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
    if ite == lambdaIte
        optPlus = lambdaPlus;
        optMinus = lambdaMinus;
    end
    
    % 全てが負だったら終了させる条件
    isOK = endOption;
    for i = 1:size(lowGap)
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
    else %更新のアルゴリズム
        lambdaMinus = [];
        lambdaPlus = [];
        UIndex = 0;
        LIndex = 0;
        if isempty(priority) == 1
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
        else %臓器ごとにpriorityがある場合
            for str = 1:strNum
                if priority(str) > 0 % 下から抑えるDVCに関して
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
                
                if priority(str) > 0 % 上から抑えるDVCに関して
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
    % stepsizeを収束に持っていくためにstepsizeを小さくする。
    stepsize = stepsize * conv;   
end
    lambda_time = toc(lambda_total);
    
    % データの書き込み
    fprintf(fileID2, 'data_Plus_x = ['); fprintf(fileID2, '%f, ', data_Plus_x); fprintf(fileID2,'];'); fprintf(fileID2,'\n'); 
    fprintf(fileID2, 'data_Minus_x = ['); fprintf(fileID2, '%f, ', data_Minus_x); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    fprintf(fileID2, 'data_Plus_y = ['); fprintf(fileID2, '%f, ', data_Plus_y); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    fprintf(fileID2, 'data_Minus_y = ['); fprintf(fileID2, '%f, ', data_Minus_y); fprintf(fileID2,'];'); fprintf(fileID2,'\n');
    
    fprintf('*** %d lambda Iteration ***\n', ite);
    lambda_time   
end
