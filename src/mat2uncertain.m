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

% single uniformly distributed randmo number in interval (0, uncertainBound)
zeta = uncertainBound * rand(1);

[~, strNum] = size(Input);
for i=1:strNum
    mat = Input{i}.mat;
    Input{i}.mat0 = mat;
    Input{i}.matPlus = mat;
    Input{i}.matMinus = mat;
    
    % the number of components
    num = numel(mat);
    numOfNoise = floor(num * rate);
    
    % select the index of matrix which is added noises
    % ind = randi(num, numOfNoise, 1); % num���Ƃ��܂��d������\��������
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
    Input{i}.matPlus = mat + uncertainBound .* abs(sp);
    
    % calculate D+ = D - uncertainBound * |D'|
    Input{i}.matMinus = mat - uncertainBound .* abs(sp);
    Input{i}.matMinus(ind(find(Input{i}.matMinus(ind)<0))) = 0;
end


end

