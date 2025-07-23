function [r, pval, U, r0_alpha] = circ_corrccnp(alpha, beta)
% [r, pval, U, r0_alpha] = circ_corrccnp(alpha, beta)
%
% Function performs a non-parametric circular-circular correlation as
% described by Zar (1999) based on Spearman rank correlation. The function
% was tested against examples given in Zar (1999) and Mardia (1975).
% Input: alpha and beta are vectors with angular data samples in radians so
%          that numel(alpha) == numel(beta)are true.
% Output: r - the correlation coefficient.
%         pval - p-value.
%         U statistic as in Mardia (1975).
%         r0_alpha - cut-off r value at alpha = 0.05 value of significance.
%
% References: Zar JH (1999). Biostatistical Analysis. 4th edition. Prentice
%               Hill.
%             Mardia, K. V. (1975). Statistics of directional data. J.
%               Royal Statist. Soc. B37: 349-393.

% By Martynas Dervinis (martynas.dervinis@gmail.com)

assert(numel(alpha) == numel(beta));

% Remove NaNs
inds = ~isnan(alpha) & ~isnan(beta);
alpha = alpha(inds);
beta = beta(inds);

% Compute ranks
n = numel(alpha);
[~, alphaSortOrder] = sort(alpha);
alphaRanks = zeros(1,n);
[~, betaSortOrder] = sort(beta);
betaRanks = zeros(1,n);
for iAB = 1:n
  alphaRanks(alphaSortOrder(iAB)) = iAB;
  betaRanks(betaSortOrder(iAB)) = iAB;
end
alphaRanks = (2*pi.*alphaRanks)./n;
betaRanks = (2*pi.*betaRanks)./n;

% Compute sum terms
sumTerm1 = 0;
sumTerm2 = 0;
sumTerm3 = 0;
sumTerm4 = 0;
for iAB = 1:n
  sumTerm1 = sumTerm1 + cos(alphaRanks(iAB) - betaRanks(iAB));
  sumTerm2 = sumTerm2 + sin(alphaRanks(iAB) - betaRanks(iAB));
  sumTerm3 = sumTerm3 + cos(alphaRanks(iAB) + betaRanks(iAB));
  sumTerm4 = sumTerm4 + sin(alphaRanks(iAB) + betaRanks(iAB));
end
sumTerm1 = sumTerm1^2;
sumTerm2 = sumTerm2^2;
sumTerm3 = sumTerm3^2;
sumTerm4 = sumTerm4^2;

% Compute the correlation coefficient
rPrime = (sumTerm1 + sumTerm2)/(n^2);
rDoublePrime = (sumTerm3 + sumTerm4)/(n^2);
r = rPrime - rDoublePrime;

% Compute the U statistic and the p-value
rVec = 0:0.000001:1;
U = 2*(n-1).*rVec;
U_pdf = exp(-0.5.*U).*(1-exp(-0.5.*U));
U_cdf = cumsum(U_pdf)/sum(U_pdf);
i = find(rVec > r, 1,'first') - 1;
pval = 1 - U_cdf(i);
U = 2*(n-1).*r;

% Compute the cut-off r value
alpha = 0.05;
r0_alpha = (-(n-1).^(-1)).*log(1 - sqrt(1 - alpha));