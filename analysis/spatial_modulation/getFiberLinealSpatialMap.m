function [fiberMaps] = getFiberLinealSpatialMap(behavior,fiber,varargin)

% USAGE
% [spatialMap] = getFiberSpatialMap(behavior,fiber,varargin)
% Calculates averaged fiber signal for a set of postions 
%
%
% Pablo Abad, NCL, 2026

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.05,@isnumeric);
addParameter(p,'nBins',100,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'linear',true);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'mode','discard'); 
addParameter(p,'order',2);
addParameter(p,'maxGap',0.1);
addParameter(p,'maxDistance',5);
addParameter(p,'sigma',2);
addParameter(p,'velocity_thresh',0.05);
addParameter(p,'min_occupancy',0);
addParameter(p,'normalize_by_occupancy',true);
addParameter(p,'plt',true);
addParameter(p,'saveFig',true);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
minTime = p.Results.minTime;
linear = p.Results.linear;
saveMat = p.Results.saveMat;
mode = p.Results.mode;
order = p.Results.order;
maxGap = p.Results.maxGap;
maxDistance = p.Results.maxDistance;
sigma = p.Results.sigma;
velocity_thresh = p.Results.velocity_thresh;
min_occupancy = p.Results.min_occupancy;
normalize_by_occupancy = p.Results.normalize_by_occupancy;
plt = p.Results.plt;
saveFig = p.Results.saveFig;

if linear
    nBins = 50;
    if isstruct(behavior)
        positions = behavior.maps;
        % maybe behavior.lin (to take all the maze and not diciding two
        % halves)
    end
else 
    nBins = nBins;
    if isstruct(behavior)
        positions{1} = [behavior.timestamps behavior.position.x behavior.position.y];
    end
end

% number of conditions
if iscell(positions)
 conditions = length(positions); 
elseif isvector(positions) || isstruct(positions)
 conditions = 1;
end
%%% TODO: conditions label
  
%% Calculate

% Interpolate fiber to position
valid_idx = fiber.timestamps >= min(behavior.timestamps) & fiber.timestamps <= max(behavior.timestamps);
fp_time_trimmed = fiber.timestamps(valid_idx);
if isfield(fiber,'red')
    fp_red_trimmed = fiber.red(valid_idx);
    fp_red_interp = interp1(fp_time_trimmed, fp_red_trimmed, behavior.timestamps, 'linear');
end
fp_green_trimmed = fiber.green(valid_idx);
fp_green_interp = interp1(fp_time_trimmed, fp_green_trimmed, behavior.timestamps, 'linear');


% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2); 

    if linear
        timestamps_position{iCond} = positions{iCond}(:,1);
        x_pos{iCond} = positions{iCond}(:,2);

        % We need to separate between left and right trials
        idx = ismember(behavior.timestamps,timestamps_position{iCond});
        fp_red{iCond} = fp_red_interp(idx);
        fp_green{iCond} = fp_green_interp(idx);


    else
        timestamps_position = positions{iCond}(:,1);
        x_pos = positions{iCond}(:,2);
        y_pos = positions{iCond}(:,3);
    end

    valid_speed_idx = v > velocity_thresh;
    
    if linear
        timestamps_position{iCond} = timestamps_position{iCond}(valid_speed_idx);
        x_pos{iCond} = x_pos{iCond}(valid_speed_idx);        
        fp_red{iCond} = fp_red{iCond}(valid_speed_idx);
        fp_green{iCond} = fp_green{iCond}(valid_speed_idx);
    else
        timestamps_position{iCond} = behavior.timestamps(valid_speed_idx);
        x_pos{iCond} = behavior.position.x(valid_speed_idx);
        y_pos{iCond} = behavior.position.y(valid_speed_idx);
        fp_red{iCond} = fp_red_interp(valid_speed_idx);
        fp_green{iCond} = fp_green_interp(valid_speed_idx);
    end
end

% Defining spatial bins
if linear
    edges = 
else
    x_edges = linspace(0, 1, nBins+1);
    y_edges = linspace(0, 1, nBins+1);
end

% Initialize maps
% if linear
% 
% else
%     map_red = zeros(nBins,nBins);
%     map_green = zeros(nBins,nBins);
%     occupancy = zeros(nBins,nBins);
%     dt_sample = median(diff(timestamps_position));
% end

% Fullfill maps
for j = 1:length(timestamps_position)
    % Initialize maps
    if linear
        map_red{j} = zeros(1,nBins);
        map_green{j} = zeros(1,nBins);
        occupancy{j} = zeros(1,nBins);
        dt_sample{j} = median(diff(timestamps_position{j}));
    else
        map_red{j} = zeros(nBins.nBins);
        map_green{j} = zeros(nBins,nBins);
        occupancy{j} = zeros(nBins,nBins);
        dt_sample{j} = median(diff(timestamps_position));
    end

    for i = 1:length(timestamps_position{j})
        x_idx = find(x_pos{j}(i) >= x_edges(1:end-1) & x_pos{j}(i) < x_edges(2:end));
        if ~linear
            y_idx = find(y_pos{j}(i) >= y_edges(1:end-1) & y_pos(i) < y_edges(2:end));
        end
    
        if ~isempty(x_idx)
            if linear
                if isfield(fiber,'red')
                    map_red{j}(x_idx) = map_red{j}(x_idx) + fp_red{j}(i);
                end
                map_green{j}( x_idx) = map_green{j}(x_idx) + fp_green{j}(i);
                occupancy{j}(x_idx) = occupancy{j}(x_idx) + dt_sample{j};

            else
                if isfield(fiber,'red')
                    map_red{j}(y_idx, x_idx) = map_red{j}(y_idx, x_idx) + fp_red{j}(i);
                end
                map_green{j}(y_idx, x_idx) = map_green{j}(y_idx, x_idx) + fp_green{j}(i);
                occupancy{j}(y_idx, x_idx) = occupancy{j}(y_idx, x_idx) + dt_sample{j};
            end
        end
    end
end

% Compute mean signal per bin
occupancy(occupancy < min_occupancy) = NaN;
if normalize_by_occupancy
    map_red_avg = map_red ./ occupancy;
    map_green_avg = map_green ./ occupancy;
else
    map_red_avg = map_red;
    map_green_avg = map_green;
end

% Smnoothing
G = fspecial('gaussian',[5 5],sigma);
if isfield(fiber,'red')
    valid_mask_red = ~isnan(map_red_avg);
    map_red_avg_zero = map_red_avg;
    map_red_avg_zero(~valid_mask_red) = 0;
    red_smoothed = imfilter(map_red_avg_zero, G, 'same');
end

valid_mask_green = ~isnan(map_green_avg);
map_green_avg_zero = map_green_avg;
map_green_avg_zero(~valid_mask_green) = 0;
green_smoothed = imfilter(map_green_avg_zero, G, 'same');




%% restructure into cell info data type

session = loadSession;

fiberMaps = [];

fiberMaps.params.smooth = sigma;
fiberMaps.params.minTime = minTime;
fiberMaps.params.nBins = nBins;
fiberMaps.params.speedThresh = speedThresh;

fiberMaps.green.occupancy = occupancy;
fiberMaps.green.count = map_green;
fiberMaps.green.z = green_smoothed;

if isfield(fiber,'red')
    fiberMaps.red.occupancy = occupancy;
    fiberMaps.red.count = map_red;
    fiberMaps.red.z = red_smoothed;
end


if saveMat
   save([session.general.name '.fiberMapAvg.mat'],'fiberMaps'); 
end

%% Plotting

if plt

    % Green
    figure;
    imagesc(fiberMaps.green.z);
    axis xy;
    colormap(jet);
    colorbar;
    title('Green fiber spatial map');
    xlabel('X bin');
    ylabel('Y bin');

    if saveFig
        saveas(gca,['SummaryFigures\fiber_green_spatial_map.png'])
    end
    
    % Red
    if isfield(fiber,'red')
        figure;
        imagesc(fiberMaps.red.z);
        axis xy;
        colormap(jet);
        colorbar;
        title('Red fiber spatial map');
        xlabel('X bin');
        ylabel('Y bin');
    end

    if saveFig
        saveas(gca,['SummaryFigures\fiber_red_spatial_map.png'])
    end
end


end