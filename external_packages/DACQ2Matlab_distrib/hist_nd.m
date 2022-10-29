function H = hist_nd(points,grid_ca)

% Bin n-d data
% points is an array of M by N points representing M observations in N dimensions
% grid_ca is a cell array of N plaid N-d arrays(from meshgrid, ndgrid) representing a grid of bin centres in each dimension
% Each element of grid_ca should be the same size
% H is the histogram of size(grid_ca{1}) in which each point is allocated to the nearest bin centre
% The delaunay triangulation is buffered between calls to prevent needless recalculation where the grid is unchanged
% However, it is the calling function's responsibility to avoid repeated calls with different grids

persistent tri lastg % To prevent recalculation of tri except where grid changes

for t = 1:numel(grid_ca)
    grid(:,t) = grid_ca{t}(:);
end

if isempty(tri) || isempty(lastg) || (numel(lastg) ~= numel(grid)) || (sum(lastg(:) ~= grid(:))>0)
    tri = delaunayn(grid,{'Qt','Qbb','Qc','Qz'});
end
lastg = grid;
k = dsearchn(grid,tri,points);
u = unique(k);
u=u(~isnan(u)); %Remove nans
H = zeros(size(grid_ca{1}));
N = hist(k,u);
if length(k) > 1 %Fix for clusters of size 1
    for t = 1:length(N)
        H(u(t)) = N(t);
    end
else
    H(u(1)) = N(N>0);
end

% ----------------------------------------------------------------------------------------------------------------------
