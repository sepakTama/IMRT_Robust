function [uncertainMat, zeta, uncertainBound] = mat2uncertain(mat, distribution)
% input: matrix and distribution to get noises
% output: uncertain matrix by reformlation some distribution
%         maximum value of zeta, that is Uncertainty
% D -> D + zeta * D'

% add noises to #(alpha * |D|) voxels
rate = 0.2;
num = numel(mat);
numOfNoise = floor(num * rate);

% select the index of matrix which is added noises
ind = randi(num, numOfNoise, 1);

% Normal Distribution
% get vector of noises
if distribution == "Normal"
    mu = 0;
    sigma = 1;
    % threshold = 1;
    vecNoise = sigma.*randn(numOfNoise,1)+mu;
    %{
    % take upper and lower bounds to noise
    vecNoise(vecNoise < -threshold) = -threshold;
    vecNoise(vecNoise >  threshold) =  threshold;
    %}
end

% Unsertainty
% calculate the bound of zeta
uncertainBound = 1;
% single uniformly distributed randmo number in interval (0, uncertainBound)
zeta = uncertainBound * rand(1);

% add noises to matrix
% D + zeta * D'
mat(ind) = mat(ind) + zeta .*vecNoise;

% change to nonnegative matrix
mat(ind(find(mat(ind)<0))) = 0;
uncertainMat = mat;

end

