function vennFromBinary(M, labels)
% plotVennFromBinary  Plot a 2-set Venn diagram from a binary matrix.
%
%   M      : [N x 2] logical or 0/1 matrix. Each column = a condition.
%   labels : 1x2 cell array with labels (optional).
%
% The function:
%   - Counts how many elements belong to each condition.
%   - Computes a proportional Venn diagram (areas ≈ counts).
%   - Computes the exact overlap area by solving for center distance.
%   - Draws the two circles and shades the overlap.
%
% MV and ChatGPT5 2025

    if nargin < 2
        labels = {'Cond 1', 'Cond 2'};
    end

    if size(M,2) ~= 2
        error('M must be an N x 2 binary matrix.');
    end

    % Count elements in each condition
    A = sum(M(:,1) ~= 0);
    B = sum(M(:,2) ~= 0);

    % Count intersection
    C = sum(M(:,1) ~= 0 & M(:,2) ~= 0);

    if C > min(A,B)
        warning('Intersection C > min(A,B); adjusting C.');
        C = min(A,B);
    end

    % Convert counts to "areas" of circles
    % Use area = pi * r^2  -> r = sqrt(area / pi)
    r1 = sqrt(A/pi);
    r2 = sqrt(B/pi);

    targetArea = C;

    % Compute distance between circle centers so that the
    % intersection area equals the number of shared elements.
    if targetArea <= 0
        % If no overlap, place circles slightly apart
        d = r1 + r2 + 0.1 * max(r1, r2);
    else
        % Function whose root corresponds to correct distance
        fun = @(d) circleIntersectionArea(d, r1, r2) - targetArea;

        % Search domain for d: 0 (full overlap) to r1+r2 (no overlap)
        d_min = 0;
        d_max = r1 + r2;

        f_min = fun(d_min);
        f_max = fun(d_max);

        if f_min * f_max > 0
            % Could not bracket the solution → fallback approximate distance
            warning('Could not bracket intersection properly. Using approximate distance.');
            d = 0.5 * (r1 + r2);
        else
            d = fzero(fun, [d_min, d_max]);
        end
    end

    % Circle centers
    center1 = [0 0];
    center2 = [d 0];

    % Draw circles
    ang = linspace(0, 2*pi, 400);
    x1 = center1(1) + r1 * cos(ang);
    y1 = center1(2) + r1 * sin(ang);
    x2 = center2(1) + r2 * cos(ang);
    y2 = center2(2) + r2 * sin(ang);

    figure; hold on
    fill(x1, y1, [0.2 0.5 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    fill(x2, y2, [1 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    axis equal off

    % Text labels
    text(center1(1), center1(2), sprintf('%s = %d', labels{1}, A), ...
         'HorizontalAlignment','center','FontSize',12);
    text(center2(1), center2(2), sprintf('%s = %d', labels{2}, B), ...
         'HorizontalAlignment','center','FontSize',12);

    % Overlap label
    text(mean([center1(1) center2(1)]), 0, sprintf('Overlap = %d', C), ...
         'HorizontalAlignment','center','FontSize',12, 'FontWeight','bold');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction: exact circle intersection area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function area = circleIntersectionArea(d, r1, r2)
% circleIntersectionArea  Exact area of intersection of two circles.
%
%   d  : distance between circle centers
%   r1 : radius of circle 1
%   r2 : radius of circle 2

    if d >= r1 + r2
        area = 0;  % no overlap
        return
    end

    if d <= abs(r1 - r2)
        % One circle completely inside the other
        area = pi * min(r1, r2)^2;
        return
    end

    % Core formula for circle-circle intersection
    alpha = 2 * acos((d^2 + r1^2 - r2^2) / (2*d*r1));
    beta  = 2 * acos((d^2 + r2^2 - r1^2) / (2*d*r2));

    area = 0.5 * r1^2 * (alpha - sin(alpha)) + ...
           0.5 * r2^2 * (beta  - sin(beta));
end
