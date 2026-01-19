function X = patchIntervals(X, t, intervals_artifacts, method)
% Patch artifacts by interpolating across specified time intervals.
%
% X: N x M  (neurons x time)
% t: 1 x M  (timestamps, monotonic)
% intervals_artifacts: K x 2  (each row [t_start t_end] in same units as t)
% method: 'linear' (default) | 'pchip' | 'nearest'
%
% Overwrites X(:, mask) with interpolated values from non-masked samples.

    if nargin < 4 || isempty(method), method = 'linear'; end
    t = t(:).';                      % 1 x M
    M = numel(t);
    N = size(X,1);

    if isempty(intervals_artifacts)
        return;
    end

    % --- Clean intervals: sort, clamp, drop invalid, merge overlaps
    iv = intervals_artifacts;
    iv = sort(iv, 2);                          % ensure start<=end per row
    iv = iv(all(isfinite(iv),2), :);           % drop NaN/Inf rows
    iv = iv(iv(:,2) >= iv(:,1), :);            % drop invalid

    if isempty(iv), return; end

    % sort by start time
    iv = sortrows(iv, 1);

    % merge overlapping/contiguous intervals
    merged = iv(1,:);
    for k = 2:size(iv,1)
        if iv(k,1) <= merged(end,2)  % overlap/contiguous (if you want a gap tolerance, adjust here)
            merged(end,2) = max(merged(end,2), iv(k,2));
        else
            merged = [merged; iv(k,:)]; %#ok<AGROW>
        end
    end
    iv = merged;

    % --- Build mask of samples to overwrite
    mask = false(1, M);
    for k = 1:size(iv,1)
        mask = mask | (t >= iv(k,1) & t <= iv(k,2));
    end

    if ~any(mask), return; end

    goodIdx = find(~mask);
    badIdx  = find(mask);

    if numel(goodIdx) < 2
        warning('Not enough non-artifact samples to interpolate.');
        return;
    end

    tg = t(goodIdx);
    tb = t(badIdx);

    % --- Patch each neuron independently
    for i = 1:N
        yi = X(i,:);
        ygood = yi(goodIdx);

        ok = ~isnan(ygood);
        if nnz(ok) < 2
            continue;
        end

        % 'extrap' avoids NaNs if intervals touch edges; change to NaN if you prefer
        ypatch = interp1(tg(ok), ygood(ok), tb, method, 'extrap');

        yi(badIdx) = ypatch;
        X(i,:) = yi;
    end
end