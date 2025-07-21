function [rho, pval, U_n] = circ_corrclnp(alpha, x)
% [rho, pval] = circ_corrclnp(alpha, x)
%
% Function performs a non-parametric circular-linear correlation as
% described by Mardia (1976). It was tested against examples in Mardia &
% Jupp (2000) and Watson (1987).
% Input: alpha - a vector sample of angular values in radians.
%        x - a vector of corresponding linear values.
% Output: rho - a non-parametric circular-linear correlation coefficient.
%         pval - a corresponding p-value.
%         U_n statistic.
%
% References: Mardia, K. V. (1976). Linear-angular correlation coefficients
%               and rhythmometry. Biometrika 63: 403-405.
%             Mardia, K. & Jupp, P. (2000). Directional Statistics. John
%               Wiley and Sons, Chichester. 
%             Watson, R. E. (1987). Two Educational Comparisons of Linear
%               and Circular Statistics. Dissertations. 2501.
%               https://ecommons.luc.edu/luc_diss/2501

% By Martynas Dervinis (martynas.dervinis@gmail.com)

% Remove NaNs
inds = ~isnan(alpha) & ~isnan(x);
alpha = alpha(inds);
x = x(inds);

% Rank values
n = numel(alpha);
[~, sortOrder] = sort(alpha);
alphaRanks = zeros(1,n);
for iA = 1:n
  alphaRanks(sortOrder(iA)) = iA;
end
alphaRanks = (2*pi.*alphaRanks)./n;

[~, sortOrder] = sort(x);
xRanks = zeros(1,n);
for iX = 1:n
  xRanks(sortOrder(iX)) = iX;
end

% Compute the D_n statistic
r = rem(n, 2);
if ~r % even
  a_n = 1/(1 + 5*((cot(pi/n))^2) + 4*((cot(pi/n))^4));
else  % odd
  a_n = (2*((sin(pi/n))^4))/((1 + cos(pi/n))^3);
end

T_c = 0; T_s = 0;
for i = 1:n
  T_c = T_c + xRanks(i)*cos(alphaRanks(i));
  T_s = T_s + xRanks(i)*sin(alphaRanks(i));
end

D_n = a_n*(T_c^2 + T_s^2);
rho = D_n;

% Compute the p-value
U_n = (24*(T_c^2 + T_s^2))/((n^2)*(n + 1));
pval = 1 - chi2cdf(U_n, 2);