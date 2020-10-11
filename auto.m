% check some result

% Parameterts
alpha = 1;
iterationNum = 5;
option = 'onlyZ';
lambdaMinus = [0;0;0;0];
lambdaPlus = [0;0;0;0];


global planC
goals = planC{planC{end}.IM}(end).IMDosimetry.goals;

% Select One of Data automatically
minID = 100;
for i=1:length(goals)
    minID = min(minID, goals(i).structNum);
end

if minID == 1
   DATA = 'Cshape';
elseif minID == 4
   DATA = 'HeadAndNeck';
elseif minID == 14
   DATA = 'Prostate';
elseif minID == 10
   DATA = 'MultiTarget';
else
   error('No Data');
end

% Select One of Data Manually
% DATA = 'Cshape'; %  1 (Core),2 (Outer Target), 
% DATA = 'HeadAndNeck'; % 4 (PTV), 5 (Cord), 6 (LtParotid), 7(RtParotid)
% DATA = 'Prostate'; %  14 (Prostate PTV), 17 (Rectum), 16 (Bladder)
% DATA = 'MultiTarget'; %  10 (Center), 11 (Inferior), 12 (Superior)
switch DATA
    case 'Cshape'
        [Input, IDs] = setTG119CShape;
    case 'HeadAndNeck'
        [Input, IDs] = setHeadAndNeck;
    case 'Prostate'
        [Input, IDs] = setProstate;        
    case 'MultiTarget'
        [Input, IDs] = setMultiTarget;
    otherwise
        error('No Data');
end

data_path = "./log/data.csv";
% data_csv = ["alpha" "iterationNum" "ratio" "t_history" "time_history"];
data_csv = [];
for ratio=0:0.01:0.4
    for UncertainBound=0.05:0.05:0.2
        for loop=1:3
            [Input, x, opt,  t_history, time_history, total_time] = UncertainIterate( ...
                            Input, alpha, iterationNum, option, lambdaMinus, lambdaPlus, ratio, UncertainBound);

            data_csv = [data_csv; [alpha iterationNum ratio UncertainBound t_history.' time_history.' total_time]];
        end
    end
end
csvwrite(data_path, data_csv);



