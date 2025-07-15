function s0 = circ_stdFisher(angularData)
% s0 = circ_stdFisher(angularData)
%   Computes the circular standard error of the mean as in Fisher & Lewis
%   (1983).
%
%   Input: angularData - a vector sample of angles in radians.
%   Output: s0.
%
%   Dependencies: circmean and circ_r of Circular Statistics Toolbox.
%
% References:
%   Fisher, N. I. & Lewis, T. (1983). Estimating the common mean direction
%   of several circular or spherical distributions with differing
%   dispersions. Biometrika 70, 333-41.

% By Martynas Dervinis (martynas.dervinis@gmail.com)

angularData = torow(angularData);
n = numel(angularData);
dataMean = circmean(angularData);
R = circ_r(angularData,[],[],2);

sumCosTerm = 0;
for i = 1:numel(angularData)
  sumCosTerm = sumCosTerm + cos(2*(angularData(i) - dataMean));
end

s0 = sqrt((1 - sumCosTerm/n)/(2*n*(R^2)));