function [LU] = IntersectIntervals(data)
% MV-ValeroLab 2024

L = data(:,1);
U = data(:,2);
% how many data intervals
n = numel(L);

[I,J] = meshgrid(1:n);

Li = L(I);
Uj = U(J);

% interval width
Wij = Uj - Li;
% ignore those with negative interval width
ignor = Wij <= 0;
Wij(ignor) = -inf;
Li(ignor) = -inf;
Uj(ignor) = inf;

% how many data intervals does each such interval lie inside?
overlapcount = sum((Li >= reshape(L,[1 1 n])) & (Uj <= reshape(U,[1 1 n])),3);

% find the sub-interval with the maximum coverage. If there is more than
% one such solution, we will choose the one with maximum interval width.
intervalscovered = max(overlapcount,[],'all');
ind = find(overlapcount == intervalscovered);

if numel(ind) == 1
    % only one interval was found with the maximum data intervals covered.
    % It must be the preferred solution.
    LU = [Li(ind),Uj(ind)];
    intervallength = Wij(ind);
else
    % choose the interval with maximum width, over those found
    [~,K] = max(Wij(ind));
    LU = [Li(ind(K)),Uj(ind(K))];
    intervallength = Wij(ind(K));
end

end