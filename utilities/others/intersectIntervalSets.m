function C = intersectIntervalSets(A, B, doMerge, tol)
%INTERSECTINTERVALSETS Intersection of two interval sets.
%
%   C = intersectIntervalSets(A, B)
%   C = intersectIntervalSets(A, B, doMerge, tol)
%
% Inputs
%   A, B     [n x 2] arrays with [start end] per row (start <= end)
%   doMerge  (optional) true to merge overlapping/contiguous results (default true)
%   tol      (optional) tolerance to treat touching intervals as contiguous (default 0)
%
% Output
%   C        [m x 2] array of intersected intervals

    if nargin < 3 || isempty(doMerge), doMerge = true; end
    if nargin < 4 || isempty(tol), tol = 0; end

    % Edge cases
    if isempty(A) || isempty(B)
        C = zeros(0,2);
        return
    end

    % Ensure proper shape
    A = double(A); B = double(B);
    if size(A,2) ~= 2 || size(B,2) ~= 2
        error('A and B must be [n x 2] interval arrays.');
    end

    % Sort by start
    A = sortrows(A, 1);
    B = sortrows(B, 1);

    % Optional: clean invalid rows
    A = A(A(:,1) <= A(:,2), :);
    B = B(B(:,1) <= B(:,2), :);

    i = 1; j = 1;
    C = zeros(0,2);

    while i <= size(A,1) && j <= size(B,1)
        a1 = A(i,1); a2 = A(i,2);
        b1 = B(j,1); b2 = B(j,2);

        % Intersection of current intervals
        s = max(a1, b1);
        e = min(a2, b2);

        if e >= s  % overlap (use e >= s for closed intervals)
            C(end+1,:) = [s e]; %#ok<AGROW>
        end

        % Move the pointer that ends first
        if a2 < b2
            i = i + 1;
        else
            j = j + 1;
        end
    end

    if doMerge && ~isempty(C)
        C = mergeIntervals(C, tol);
    end
end

function M = mergeIntervals(X, tol)
%MERGEINTERVALS Merge overlapping or contiguous intervals in X.
%   tol allows merging if next.start <= current.end + tol

    X = sortrows(X,1);
    M = X(1,:);
    for k = 2:size(X,1)
        if X(k,1) <= M(end,2) + tol
            M(end,2) = max(M(end,2), X(k,2));
        else
            M(end+1,:) = X(k,:); %#ok<AGROW>
        end
    end
end
