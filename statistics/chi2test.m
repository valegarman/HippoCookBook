
function [p, stats] = chi2test(varargin)    
% chi2test - Compute chi2 test
%
% USAGE
%    [p, stats] = chi2test(varargin)
%    [p, stats] = chi2test([obs1, total1],[obs2, total2],etc)
% 
% INPUTS
% x1, x2... xn, being xi a 2 elements column vector on the form [observation,
%   total]. Example chi2stest([41, 141],[10,50]) for a total of 141 and 50
%   observations, respectively
% OUPTUT
%
% P
% stats.
% 
% By Manu Valero 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:nargin
        n(ii) = varargin{ii}(1);
        N(ii) = varargin{ii}(2);
    end
    k = length(n);
    
    % Pooled estimate of proportion
    p0 = sum(n)/sum(N);
    
    % Expected counts under H0 (null hypothesis)
    n0 = N * p0;
    
    % Chi-square test
    observed = [n(:)' (N-n)];
    expected = [n0(:)' (N-n0)];
    
   [~,p,stats] = chi2gof([1:length(observed)],'freq',observed,'expected',expected,'ctrs',[1:length(observed)],'nparams',k);
   stats.p = p;
   stats.testName = 'χ2 Goodness of Fit Test'; 
   
   tbl = [];
   if nargin > 2   
        comb=nchoosek(1:nargin,2);
        for ii = 1:size(comb,1)
            idx = comb(ii,:);
            observed = [n(idx) (N(idx)-n(idx))];
            expected = [n0(idx) (N(idx)-n0(idx))];
            [~,pp,~] = chi2gof([1:length(observed)],'freq',observed,'expected',expected,'ctrs',[1:length(observed)],'nparams',2);
            
            tbl = [tbl; comb(ii,:) pp];
            fprintf('%3.i %3.i      %3.3e\n',comb(ii,1), comb(ii,2), pp); %\n
        end
        stats.tbl = tbl;
   end
   
   
end