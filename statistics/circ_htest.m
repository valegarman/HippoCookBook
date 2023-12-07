function [pval F] = circ_htest(alpha1,alpha2,r1,r2)

%   Parametric Hotelling paired sample test for equal angular means. Can be 
%   used as a paired-sample test for circular data. The test has 2 by N-2
%   degrees of freedom where N is the number of pairs in the phase 
%   distribution. 
%
%   If the individual phase values in alpha are mean values of 
%   distributions, this concearns a second-level analysis. In that case,
%   the associated average vector lengths of these distributions are taken 
%   into account.
%
%   H0: the two populations have equal mean phase angles
%   HA: the two populations have unequal mean phase angles
%   
%   Usage:
%     [pval, F] = circ_htest(alpha1, alpha2)
%     [pval, F] = circ_htest(alpha1, alpha2, r1 ,r2)
%
%   Input:
%     alpha1:  angles in radians or degrees
%     alpha2:  angles in radians or degrees
%     r1: (optional) resultant vector lengths of alpha1, in case the values
%     in alpha1 consist of mean values of distributions themselves (i.e. 
%     second-level analysis)
%     r2: (optional) resultant vector lengths of alpha2, in case the values
%     in alpha2 consist of mean values of distributions themselves (i.e. 
%     second-level analysis)
%
%   Notes:
%     Inputs can be either a row or column vector. Inputs alpha1 and alpha2 
%     need to have the same length. Optional r1 and r2 need to be the same 
%     length as alpha, and each value in r should have the same index as 
%     their alpha counterpart.
%
%   Output:
%     pval    p-value of the Hotelling paired-sample test. Discard H0 if
%             pval is small.
%     F       the test statistic of the Hottelling test
%
% References:
%   Biostatistical Analysis, J. H. Zar (1999)
%
% RL van den Brink, 2014
% r.l.van.den.brink@fsw.leidenuniv.nl

%% check the input
if nargin == 0
    error('not enough input arguments')
elseif nargin == 1
    error('provide at least two phase distributions')
elseif nargin == 2
    meanofmeans = 0; %this isn't a second-level analysis
    if length(alpha1) ~= length(alpha2)
        error('phase distributions must be the same length')
    end
elseif nargin == 3
    error('provide at least to vecotrs of mean resultant length of the phase distributions')
elseif nargin == 4
    meanofmeans = 1; %his is a second level analysis
    if length(r1) ~= length(alpha1) || length(r2) ~= length(alpha2) || length(r1) ~= length(alpha2) || length(r2) ~= length(alpha1)
        error('vector of resultant lengths must be the same length as phase values')
    end    
    %make sure they are row vectors
    if size(r1,2) > size(r1,1); r1 = r1'; end
    if size(r2,2) > size(r2,1); r2 = r2'; end
else
    error('too many input arguments')
end

%% transform and normalize if needed
%make sure these are all row vectors
if size(alpha1,2) > size(alpha1,1); alpha1 = alpha1'; end
if size(alpha2,2) > size(alpha2,1); alpha2 = alpha2'; end

%test if phases are radians or degrees, and normalize to between 0 and 2pi
%rad if needed (for the calculations below the phases have to be positive
%and in radians)
alpha = alpha1(~isnan(alpha1)); %ignore NaN
m = min(alpha(:));
M = max(alpha(:));

if m >=0 && M <= 360,
	alpha1 = ang2rad(alpha1);
elseif m >= -pi && M <= pi
	alpha1(alpha1 < 0) = alpha1(alpha1 < 0) + (2*pi);
elseif m >= -180 && M <= 180,
    alpha1(alpha1 < 0) = alpha1(alpha1 < 0) + 360;
    alpha1 = ang2rad(alpha1);
end

alpha = alpha2(~isnan(alpha2)); %ignore NaN
m = min(alpha(:));
M = max(alpha(:));

if m >=0 && M <= 360,
	alpha2 = ang2rad(alpha2);
elseif m >= -pi && M <= pi
	alpha2(alpha2 < 0) = alpha2(alpha2 < 0) + (2*pi);
elseif m >= -180 && M <= 180,
    alpha2(alpha2 < 0) = alpha2(alpha2 < 0) + 360;
    alpha2 = ang2rad(alpha2);
end

clear alpha

%% Run the test
k = length(alpha1); %number of phase pairs

% take into account the resultant vector lengths if this is a second-level
% analysis (i.e. the this is a test of the mean of means)
if meanofmeans
    %get normalized rectangular coordinates
    X = r2.*cos(alpha2) - r1.*cos(alpha1);
    Y = r2.*sin(alpha2) - r1.*sin(alpha1);
else
    %get rectangular coordinates
    X = cos(alpha2) - cos(alpha1);
    Y = sin(alpha2) - sin(alpha1);
end

%get (normalized) grand mean rectangular coordinates
Xbar = mean(X);
Ybar = mean(Y);

sigmaxsq = sum(X.^2) - (sum(X)^2 ./ k);
sigmaysq = sum(Y.^2) - (sum(Y)^2 ./ k);
sigmaxy  = sum(X.*Y) - ((sum(X).*sum(Y)) / k);

%get F-statistic, see Biostatistical Analysis, Zar (1999), section 27.13 for equation
F = ((k*(k-2)) / 2) *  (((Xbar^2*sigmaysq) - (2*Xbar * Ybar * sigmaxy) + (Ybar^2 * sigmaxsq)) / ...
                        ((sigmaxsq * sigmaysq) - (sigmaxy^2)) );
%get p-value                    
pval = 1 - fcdf(F,2,k-2); 

%converts degrees to radians
function alpha = ang2rad(alpha)
alpha = alpha * pi /180;

% %converts radians to degrees
% function alpha = rad2ang(alpha)
% alpha = alpha / pi *180;
