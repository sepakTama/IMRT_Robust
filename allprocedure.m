% all procedure

% Parameterts
alpha = 1;
iterationNum = 5;
lambdaIte = 1;
option = 'onlyZ';
priority = [1;0];

% alphaÇÃëùâ¡ó 
delta = 0.01;
end_alpha = 1;

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

global init;
init = [];
global init_lambdaMinus;
global init_lambdaPlus;


colNum = floor((end_alpha - alpha)/delta);
rowNum = 2;

[~, beamletNum] = size(Input{1}.mat);
[~, strNum] = size(Input);
for str = 1:strNum
    if isempty(Input{str}.up) ~= 1
        rowNum = rowNum + max(size(Input{str}.up));
    end
    if isempty(Input{str}.lp) ~= 1
        rowNum = rowNum + max(size(Input{str}.lp));
    end
end
rowNum = rowNum + beamletNum;
data_csv = zeros(rowNum, colNum);

iteCol = 1;
    
for mu = alpha:delta:end_alpha
    [Input, x, optPlus, optMinus, lambda_time] = randomIteration(Input, mu, iterationNum, lambdaIte, option, priority);
    
    data_csv(:, iteCol) = [mu; optPlus; optMinus; lambda_time; x];
    iteCol = iteCol + 1;
end
csvwrite('HeadAndNeck.dat', data_csv);
init = [];
