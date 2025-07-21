
function [p, z] = getFisherZtransformationTest(A, B, an, bn)

% MV-NCL 2025

if any(size(A) ~= size(B))
    error('Both coef matrix inputs should have the same size...');
end

if length(an) == 1
    an = an .*ones(size(A));
end

if length(bn) == 1
    bn = bn .*ones(size(B));
end

p = nan(size(A));
z = nan(size(A));
for ii = 1:size(A,1)
    for jj = 1:size(A,2)
        r1 = A(ii,jj);    % correlation from sample 1
        n1 = an(ii,jj);   % sample size 1
        r2 = B(ii,jj);    % correlation from sample 2
        n2 = bn(ii,jj);   % sample size 2
        
        % Fisher's r-to-z transformation
        z1 = 0.5 * log((1 + r1) / (1 - r1));
        z2 = 0.5 * log((1 + r2) / (1 - r2));

        % Standard error
        SE = sqrt(1 / (n1 - 3) + 1 / (n2 - 3));

        % z statistic
        z_value = (z1 - z2) / SE;

        % Two-tailed p-value
        p_value = 2 * (1 - normcdf(abs(z_value)));

        p(ii,jj) = p_value;
        z(ii,jj) = z_value;
    end
end


end