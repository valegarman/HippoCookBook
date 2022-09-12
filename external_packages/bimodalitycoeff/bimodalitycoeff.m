% +------------------------------------------------------+
% |          Bimodality Coefficient Calculation          |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        12/28/20 | 
% +------------------------------------------------------+
% 
% function: [BF, BC] = bimodalitycoeff(x)
%
% Input:
% x - data in the time domain; x could be vector or matrix
%     with time across columns and indexes across rows
% 
% Output:
% BF - bimodality flag (treshold of 5/9)
% BC - bimodality coefficient (unbiased)

function [BF, BC] = bimodalitycoeff(x)

% check if x is vector and if it is - 
% represent it as a column-vector
if isvector(x), x = x(:); end

% determine the data size along its first dimension
N = size(x, 1);

% determine the data skewness (unbiased)
S = skewness(x, 0);

% determine the data kurtosis (unbiased)
% (the normal distribution has kurtosis of zero)
K = kurtosis(x, 0) - 3;

% calculate the bimodality coefficient (unbiased)
BC = (S.^2 + 1) ./ (K + 3*(N-1)^2/(N-2)/(N-3));

% determine the bimodality flag (using +5% margin)
BF = BC > 1.05*5/9;

end