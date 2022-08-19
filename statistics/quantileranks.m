function Y = quantileranks(X, q, permtie)
% return quantile ranks of the values in X based on
% their sorted order.
%
% X is a vector of values, possibly with NaNs.
% q determines quantile class i.e.
%   q =   4, quartiles   ( 4-quantile)
%   q =  10, deciles     (10-quantile)
%   q = 100, percentiles (100-quantile)
%   ...
% permties is a boolean flag to shuffle ranks over
% duplicate values, ties.
% Y is a vector of quantile ranks for the corresponding X values.
%
% notes:
% 1. as possible, the values in X are distributed evenly over quantile
% rank groups, could be verified using tabulate(Y).
% 2. NaNs are assigned the rank of zero.
% 3. duplicate values, ties, are distributed randomly over
% different groups (if assigned different ranks).
% 4. to find all values in the ith quantile group, e.x. the 5th decile,
% use X(Y==grp#), e.x. X(Y==5).
%
% Author: Mohammad Rawashdeh
% rawashmy at mail dot uc dot edu
% April 2013.
if nargin == 2
    permtie = true;
elseif nargin ~= 3
    msg = 'ERROR: the function at least takes two inputs.';
    error('quantileranks:INPUT_ERROR', msg);
end

Y = zeros(size(X));
if ~isvector(X)
    msg = 'WARNING: first input must be a vector; returning ranks of zeros.';
    warning('quantileranks:INPUT_ERROR', msg);
    return;
end

non_nans_n = sum(~isnan(X)); % count of non-nans

% quantile rank positions in sorted X
q_pos = zeros(q,1);
for rank = 1:q
    q_pos(rank) = ceil(rank/q * non_nans_n);
end

[~,IX] = sort(X);

% initialize
rank = 1;
IX_first = 1;
IX_last = q_pos(rank);

% assign ranks
while (rank <= q)
    % pick highest possible rank (in case q > non_nans_n)
    while ( (rank < q) && (q_pos(rank+1) == IX_last) )
        rank = rank+1;
    end
    % assign ranks
    for k = IX_first:IX_last
        Y(IX(k)) = rank;
    end
    % next rank
    IX_first = IX_last + 1;
    rank = rank + 1;
    if (rank <= q)
        IX_last = q_pos(rank);
    end
end

% break ties
if permtie
    YT = Y(IX);
    k = 1;
    rank = YT(k);
    while (k <= non_nans_n)
        % find group boundary
        while (k < non_nans_n && YT(k+1) == rank)
            k = k + 1;
        end
        
        % if there is a tie at boundary then break
        % the tie by shuffling the ranks over tie values.
        if ( (k < non_nans_n) && (X(IX(k)) == X(IX(k+1))) )
            tie_value = X(IX(k));
            % find tie first occurrence, move backward
            bkw = k;
            while (bkw > 1 && X(IX(bkw-1)) == tie_value)
                bkw = bkw - 1;
            end
            
            % find tie last occurrence, move forward
            frw = k + 1;
            while (frw < non_nans_n && X(IX(frw+1)) == tie_value)
                frw = frw + 1;
            end
            
            % shuffle and reassign ranks
            IX_tie = IX(bkw:frw);
            shuffle = randperm(length(IX_tie));
            Y(IX_tie) = Y(IX_tie(shuffle));
            
            k = frw + 1;
        else
            k = k + 1;
        end
        
        if (k <= non_nans_n)
            rank = YT(k);
        end
    end
end

end