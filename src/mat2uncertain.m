function [Input, uncertainBound] = mat2uncertain(Input, distribution, ratio, uncertainBound)
% input: matrix and distribution to get noises
% output: uncertain matrix by reformlation some distribution
%         maximum value of zeta, that is Uncertainty

% add noises to #(alpha * |D|) voxels
%rate = 0.2;
rate = ratio;

% Unsertainty
% calculate the bound of zeta
%uncertainBound = 0.1;

global planC
goals = planC{planC{end}.IM}(end).IMDosimetry.goals;
% Select One of Data automatically
minID = 100;
for i=1:length(goals)
    minID = min(minID, goals(i).structNum);
end
if minID == 1
    DATA = 'CShape';
elseif minID == 4
    DATA = 'HeadAndNeck';
elseif minID == 10
    DATA = 'MultiTarget';
elseif minID == 14
    DATA = 'Prostate';
end

% 1 -> PTV; 0 -> OAR
switch DATA
    case 'CShape'
        ptv = [1; 0];
    case 'HeadAndNeck'
        ptv = [1; 0; 0; 0];
    case 'MultiTarget'
        ptv =  [1; 1; 1];
    case 'Prostate'
        ptv = [1; 0; 0];
end
% boolean to use fix param or not, 1 -> true 0 -> bool
bool_use_fix = 0;
% Fix noises in PTV or OAR to check the movement of one of both


[~, strNum] = size(Input);
for i=1:strNum
    mat = Input{i}.mat;
    Input{i}.mat0 = mat;
    Input{i}.matPlus = mat;
    Input{i}.matMinus = mat;
    num = numel(mat);
    
    if ~ptv(i) && bool_use_fix % the condition that structure is OAR and bool is true
        % fix parameter
        fix_rate = 0.1;
        bound = 0.1;
        % the number of components
        numOfNoise = floor(num * fix_rate);
    else % the condition that structure is OAR or bool is false
        % the number of components
        numOfNoise = floor(num * rate);
        bound = uncertainBound;
    end
    
    % single uniformly distributed randmo number in interval (0, uncertainBound)
    zeta = bound * rand(1);
    
    % select the index of matrix which is added noises
    % ind = randi(num, numOfNoise, 1); % num?½?½?½Æ‚ï¿½?½Ü‚ï¿½?½d?½?½?½?½?½?½Â”\?½?½?½?½?½?½?½?½
    ind = transpose(randperm(num, numOfNoise));
    ind = sort(ind);
    
    % Normal Distribution
    % get vector of noises
    if distribution == "Normal"
        mu = 0;
        sigma = 1;
        % lambda function
        normal = @(numOfNoise, mu, sigma) sigma.*randn(numOfNoise,1)+mu;
        vecNoise = normal(numOfNoise, mu, sigma);
    end
    % make sparse matrix adding noise
    [r_size, c_size] = size(mat);
    sp = sparse(r_size, c_size);
    sp(ind) = vecNoise;
    
    % calculate D0 = D + zeta * D'
    Input{i}.mat0 = mat + zeta .* sp;
    Input{i}.mat0(ind(find(Input{i}.mat0(ind)<0))) = 0;

    vecNoise = abs(vecNoise);
    % calculate D+ = D + uncertainBound * |D'|
    % no need to consicder non-negative constraints
    Input{i}.matPlus = mat + bound .* abs(sp);

    % calculate D+ = D - uncertainBound * |D'|
    Input{i}.matMinus = mat - bound .* abs(sp);
    Input{i}.matMinus(ind(find(Input{i}.matMinus(ind)<0))) = 0;
end


end

