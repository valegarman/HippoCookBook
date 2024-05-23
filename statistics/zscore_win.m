function [z,sigma,mu] = zscore_win(x,binarized_win)
% zscore for matrices and vectors based on the std and the mean estimated
% from a given interval

% INPUTS
% X                 M x N matrix 
% binarized_win     Binarized vector indicating as true (or 1) the target values for computing mu and sigma
%
% NCL-MV 2024

if nargin < 2
    binarized_win = ones(size(x,2),1);
end

if size(x,1) ~= length(binarized_win)
    x = x';
end

if size(x,1) ~= length(binarized_win)
    error('Input data and binarized window has not the same number of elements!');
end

z = zeros(size(x));
for ii = 1:size(x,2)
    z(:,ii) = (x(:,ii) - mean(x(binarized_win,ii)))./std(x(binarized_win,ii));
    mu(ii) = mean(x(binarized_win,ii));
    sigma(ii) = std(x(binarized_win,ii));
end

end
