% check information of matrix
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

[~, strNum] = size(Input);
for str = 1:strNum
   fprintf('===Structure{%d}===', str);
   num = numel(Input{str}.mat);
   %med = median(Input{str}.mat);
   me = mean(Input{1}.mat, [1, 2]);
   nnz(Input{str}.mat);
   fprintf('number of elements = {%d}', num);
   %fprintf('median = {%f}', med);
   %fprintf('mean = {%f}', me);
   %fprintf('number of nonzero = {%d}', nn);
end
