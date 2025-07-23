function [t0, t] = circ_confmeanFisher(alpha, xi, dim)
% [t0, t] = circ_confmeanFisher(alpha, xi, dim)
%   Computes confidence limits on the mean for circular data without
%   making assumptions about the data distribution.
%
%   Input: alpha - sample of angles in radians. If it is a matrix, then
%                  confidence intervals are calculated along the specified
%                  dimension.
%          xi - (1-xi)-confidence limits are computed, default 0.05.
%          dim  - compute along this dimension, default is 1.
%
%   Output: t0 - half (1-xi)% confidence interval.
%           t - mean +- d yields upper/lower (1-xi)% confidence limit.
%
%   Dependencies: circ_std of Circular Statistics Toolbox.
%
% References:
%   Fisher, N. I. & Lewis, T. (1983). Estimating the common mean direction
%   of several circular or spherical distributions with differing
%   dispersions. Biometrika 70, 333-41.

% By Martynas Dervinis (martynas.dervinis@gmail.com)

if nargin < 3
  dim = 1;
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
  xi = 0.05;
end

if isempty(alpha)
  t0 = NaN;
  t = [NaN; NaN];
  return
end

% compute the mean and s0
alphaSize = size(alpha);
if numel(alphaSize) > 2
  error('Only vectors and 2-D matrices are accepted as input.');
end
if alphaSize(1) > 1 && alphaSize(2) > 1
  if dim == 1
    oppositeDim = 2;
  elseif dim == 2
    oppositeDim = 1;
  else
    error('Can only calculate confidence intervals along dimension 1 or 2.');
  end
  nCond = size(alpha,oppositeDim);
  u = nan(1, nCond);
  s0 = nan(1, nCond);
  for iCond = 1:nCond
    if dim == 1
      u(iCond) = circmean(alpha(~isnan(alpha(:,iCond)),iCond));
      s0(iCond) = circ_std(alpha(~isnan(alpha(:,iCond)),iCond));
    elseif dim == 2
      u(iCond) = circmean(alpha(iCond,~isnan(alpha(iCond,:))));
      s0(iCond) = circ_stdFisher(alpha(iCond,~isnan(alpha(iCond,:))));
    end
  end
else
  u = circmean(alpha(~isnan(alpha)));
  s0 = circ_stdFisher(alpha(~isnan(alpha)));
end

% Compute u_alpha
u_alpha = norminv(1-xi/2);

% Compute confidence intervals
t0 = asin(u_alpha.*s0);
for iT = 1:numel(t0)
  if ~isreal(t0(iT))
    t0(iT) = NaN;
  end
end
t = [u + t0; u - t0];