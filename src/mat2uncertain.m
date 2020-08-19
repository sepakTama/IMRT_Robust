function [CShape, uncertainBound] = mat2uncertain(CShape, distribution)
% input: matrix and distribution to get noises
% output: uncertain matrix by reformlation some distribution
%         maximum value of zeta, that is Uncertainty

% add noises to #(alpha * |D|) voxels
rate = 0.2;

% Unsertainty
% calculate the bound of zeta
uncertainBound = 0.1;

% single uniformly distributed randmo number in interval (0, uncertainBound)
zeta = uncertainBound * rand(1);

[~, strNum] = size(Input);
for i=1:strNum
    mat = CShape{i}.mat;
    CShape{i}.mat0 = mat;
    CShape{i}.matPlus = mat;
    CShape{i}.matMinus = mat;
    
    % the number of components
    num = numel(mat);
    numOfNoise = floor(num * rate);
    
    % select the index of matrix which is added noises
    % ind = randi(num, numOfNoise, 1); % numÇæÇ∆Ç§Ç‹Ç≠èdï°Ç∑ÇÈâ¬î\ê´Ç™Ç†ÇÈ
    ind = transpose(randperm(num, numOfNoise));
    
    % Normal Distribution
    % get vector of noises
    if distribution == "Normal"
        mu = 0;
        sigma = 1;
        % lambda function
        normal = @(numOfNoise, mu, sigma) sigma.*randn(numOfNoise,1)+mu;
        vecNoise = normal(numOfNoise, mu, sigma);
    end
    
    % calculate D0 = D + zeta * D'
    CShape{i}.mat0(ind) = mat(ind) + zeta .* vecNoise;
    CShape{i}.mat0(ind(find(CShape{i}.mat0(ind)<0))) = 0;
    
    vecNoise = abs(vecNoise);
    % calculate D+ = D + uncertainBound * |D'|
    % no need to consicder non-negative constraints
    CShape{i}.matPlus(ind) = mat(ind) + uncertainBound .* vecNoise;
    
    % calculate D+ = D - uncertainBound * |D'|
    CShape{i}.matMinus(ind) = mat(ind) - uncertainBound .* vecNoise;
    CShape{i}.matMinus(ind(find(CShape{i}.matMinus(ind)<0))) = 0;
end


end

