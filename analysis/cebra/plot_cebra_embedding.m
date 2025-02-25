function [outputArg1,outputArg2] = plot_cebra_embedding(embedding,label,varargin)

p = inputParser;

addParameter(p,'gray',false);
addParameter(p,'idx_order',[1 2 3]);
addParameter(p,'pointSize',1);

parse(p,varargin{:});

gray = p.Results.gray;
idx_order = p.Results.idx_order;
pointSize = p.Results.pointSize;

% define right and left indices

r_ind = label(:,2) == 1;
l_ind = label(:,3) == 1;

cmap = colormap(flip(brewermap(100,'RdYlBu')));

% Set colors based on 'gray' parameter
if ~gray
    r_cmap = cmap; % MATLAB equivalent to 'rainbow'
    l_cmap = cmap;
    r_c = label(r_ind, 1);
    l_c = label(l_ind, 1);
else
    r_cmap = gray;
    l_cmap = gray;
    r_c = [0.5, 0.5, 0.5]; % Use gray color
    l_c = [0.5, 0.5, 0.5];
end

% Define indices order
idx1 = idx_order(1);
idx2 = idx_order(2);
idx3 = idx_order(3);

% Plot right side
scatter3(embedding(r_ind, idx1), ...
             embedding(r_ind, idx2), ...
             embedding(r_ind, idx3), ...
             pointSize, r_c, 'filled');
colormap(r_cmap);

% Plot left side
hold on;
scatter3(embedding(l_ind, idx1), ...
             embedding(l_ind, idx2), ...
             embedding(l_ind, idx3), ...
             pointSize, l_c, 'filled');
colormap(l_cmap);

% Set axis properties
grid('off');
axesHandles = findall(gcf, 'Type', 'axes');
ax = axesHandles(1);
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';

% Make panes transparent
ax.BoxStyle = 'full';
ax.XColor = 'none';
ax.YColor = 'none';
ax.ZColor = 'none';

hold off;


end