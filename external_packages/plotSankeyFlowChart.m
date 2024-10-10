function plotSankeyFlowChart(data,options)
%% Load default options if not provided
try
    color_map = options.color_map;
catch
    color_map = 'parula';
end
try
    flow_transparency = options.flow_transparency;
catch
    flow_transparency = 0.2;
end

try
    bar_width = options.bar_width;
catch
    bar_width = 20;
end

try
    show_perc = options.show_perc;
catch
    show_perc = true;
end

try
    text_color = options.text_color;
catch
    text_color = 'w';
end

try
    show_cat_labels = options.show_cat_labels;
catch
    show_cat_labels = true;
end

try
    show_layer_labels = options.show_layer_labels;
catch
    show_layer_labels = true;
end

try
    show_legend = options.show_legend;
catch
    show_legend = true;
end


if istable(data)
    input_data = nan(height(data),width(data));
    for t =1:width(data)
        col = table2array(data(:, t));
        if iscategorical(col) || ischar(col)
            [~, ia, col] = unique(col,'stable');
            layer_cat_names{t} = cellstr(table2array(data(ia,t)));
        end
        input_data(:,t) = col;
    end
else
    for t =1:size(data,2)
        layer_cat_names{t} = cellstr(num2str(unique(data(:,t))));
    end
    input_data = data;
end

n_obs = size(input_data,1);
n_layers = size(input_data,2);

%% Detect categories used:
categories = sort(unique(input_data),'ascend');
max_n_cats_layer = length(categories);
cat_names = [];
if ischar(color_map)
    all_colors = eval([color_map '(n_layers * max_n_cats_layer)']);
else
    if sum([n_layers * max_n_cats_layer 3] == size(color_map)) == 2
        all_colors = color_map;
    else
        warning('input color map shoud have n_layers * max_n_cats_layer rows and 3 columns, using parula instead');
        all_colors = parula(n_layers * max_n_cats_layer);
    end
end

%% stacked bars for each category at each layer (t):
tidx = 0;
cidx = 1;
for t = 1:n_layers
    if istable(data)
        layer_point_names{t} = data.Properties.VariableNames{t};
    else
        layer_point_names{t} = ['Layer ',num2str(t)];
    end
    barcolors{t} = all_colors(tidx+1:tidx+max_n_cats_layer,:);
    tidx = tidx + max_n_cats_layer;
    cidx_t = 1;
    for c = 1:max_n_cats_layer
        bars{t}(1,c) = sum(input_data(:,t)==categories(c)); % sum categories...
        if bars{t}(1,c) > 0
            if istable(data)
                cat_names{cidx} = layer_cat_names{t}{cidx_t};
            else
                cat_names{cidx} = ['cat' num2str(cidx_t) '-layer' num2str(t)];
            end
            cidx_t = cidx_t + 1;
        else
            cat_names{cidx} = '-';
        end
        cidx = cidx + 1;
    end
end


%% change between layer t and t+1 for each category:
for t = 1:n_layers-1
    for c1 = 1:max_n_cats_layer % change from 1 category at t...
        for c2 = 1:max_n_cats_layer % to all other categories at t+1
            change{t}(c1,c2) = sum(input_data(:,t) == categories(c1) & input_data(:,t+1) == categories(c2));
        end
    end
end


%% Create the Figure
% layers on horizontal axis:
X = 0:n_layers-1;
figure();
y1_category_points=[];

for t=1:n_layers-1
    y1_category_points = sankey_alluvialflow(bars{t}, bars{t+1}, change{t}, X(t), X(t+1), y1_category_points,barcolors{t},barcolors{t+1},flow_transparency,bar_width);
end


%% Add text (% in each category and layer)
ymax = 1.1 * n_obs;
gap = (ymax-n_obs) / (max_n_cats_layer-1);

for t = 1:n_layers
    for c = 1:max_n_cats_layer
        val = roundn(bars{t}(c)/n_obs*100,0);
        if show_perc
            if val > 0
                if show_cat_labels
                    txt2print = [layer_cat_names{t}{c} ' [' num2str(val) '%]'];
                else
                    txt2print = [num2str(val) '%'];
                end
                text(X(t), bars{t}(c)/2 + sum(bars{t}(1:c-1)) + (c-1)*gap, txt2print,'HorizontalAlignment','center','Color',text_color)
            end
        else
            if show_cat_labels
                if val > 0
                    text(X(t), bars{t}(c)/2 + sum(bars{t}(1:c-1)) + (c-1)*gap, layer_cat_names{t}{c},'HorizontalAlignment','center','Color',text_color)
                end
            end
        end
    end
    if show_layer_labels
        text(X(t), bars{t}(c) + sum(bars{t}(1:c-1)) + (c+1)*gap, layer_point_names{t},'HorizontalAlignment','center','Interpreter','None')
    end
end

%% Legend
if show_legend
    final_legend = [];
    lidx = 1;
    for c = 1:length(cat_names)
        if ~strcmp(cat_names{c},'-')
            h(lidx) = scatter(nan, 0, 's', 'MarkerFaceColor', all_colors(c,:),'MarkerEdgeColor', 'k');
            final_legend{lidx} = cat_names{c};
            lidx = lidx + 1;
        end
    end
    legend(h, final_legend,'Interpreter','none','location','eastoutside');
end
end

function h = sankey_alluvialflow(Bars1, Bars2, change, x1, x2,last_category_points,BarColors1,BarColors2,ChangeTransparancy,BarWidth)


Height = sum(Bars1)*1.1;
Gap = sum(Bars1)*0.1 / (length(Bars1)-1);
axis ij % origin is top left
axis off

hold on

% These are the top points for each left category, with gaps added.
if isempty(last_category_points)
    y1_category_points = [0 cumsum(Bars1)] + (0:numel(Bars1)) .* Gap;
    y1_category_points(end) = [];
else
    y1_category_points=last_category_points;
end

% These are the top points for each right category, with gaps added.
y2_category_points = [0 cumsum(Bars2)] + (0:numel(Bars2)) .* Gap;
y2_category_points(end) = [];
h=y2_category_points;


% Draw the patches, an entire left category at a time
right_columns_so_far = y2_category_points(1:end); % Start at the beginning of each right category and stack as we go.
patches_per_left_category = size(change, 2);
for k_left = 1:size(change, 1) % for each row
    
    % Calculate the coordinates for all the patches split by the
    % Split the left category
    left_patch_points = [0 cumsum(change(k_left, :))] + y1_category_points(k_left);
    patch_top_lefts = left_patch_points(1:end-1);
    patch_bottom_lefts = left_patch_points(2:end);
    
    % Compute and stack up slice of each right category
    patch_top_rights = right_columns_so_far;
    patch_bottom_rights = patch_top_rights + change(k_left, :);
    right_columns_so_far = patch_bottom_rights;
    
    % Plot the patches
    
    % X coordinates of patch corners
    [bottom_curves_x, bottom_curves_y] = get_curves(x1+0.1, patch_bottom_lefts, x2-0.1, patch_bottom_rights);
    [top_curves_x,    top_curves_y]    = get_curves(x2-0.1, patch_top_rights,   x1+0.1, patch_top_lefts);
    X = [ ...
        repmat([x1; x1], 1, patches_per_left_category); % Top left, bottom left
        bottom_curves_x;
        repmat([x2; x2], 1, patches_per_left_category); % Bottom right, top right
        top_curves_x
        ];
    
    % Y coordinates of patch corners
    Y = [ ...
        patch_top_lefts;
        patch_bottom_lefts;
        bottom_curves_y;
        patch_bottom_rights;
        patch_top_rights;
        top_curves_y
        ];
    
    patch('XData', X, 'YData', Y, 'FaceColor', BarColors1(k_left,:), 'FaceAlpha', ChangeTransparancy, 'EdgeColor', 'none');
end % for each row

% plot left category bars
for i=1:numel(y1_category_points)
    y1=[y1_category_points; (y1_category_points + Bars1)];
    plot(ones(2, 1)*x1, y1(:,i), 'Color', BarColors1(i,:),'LineWidth',BarWidth);
end
hold on

% plot right category bars
for i=1:numel(y2_category_points)
    y2=[y2_category_points; (y2_category_points + Bars2)];
    plot(ones(2, 1)*x2, y2(:,i), 'Color', BarColors2(i,:),'LineWidth',BarWidth);
end


end % alluvialflow

function [x, y] = get_curves(x1, y1, x2, y2)
% x1, x2: scalar x coordinates of line start, end
% y1, y2: vectors of y coordinates of line start/ends
Npoints = 15;
t = linspace(0, pi, Npoints);
c = (1-cos(t))./2; % Normalized curve

Ncurves = numel(y1);
y = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints,1) .* repmat(c', 1, Ncurves);
x = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end  % get_curve
