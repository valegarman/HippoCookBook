function [fiberMaps] = getFiberSpatialMap(behavior,fiber,varargin)

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
addParameter(p,'nBins',75,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'mode','discard'); 
addParameter(p,'order',2);
addParameter(p,'maxGap',0.1);
addParameter(p,'maxDistance',5);
addParameter(p,'sigma',2);
addParameter(p,'velocity_thresh',0);
addParameter(p,'min_occupancy',0);
addParameter(p,'normalize_by_occupancy',true);
addParameter(p,'plt',true);
addParameter(p,'saveFig',true);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
minTime = p.Results.minTime;
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


if isstruct(behavior)
    positions{1} = [behavior.timestamps behavior.position.x behavior.position.y];
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

    valid_speed_idx = v > velocity_thresh;

    timestamps_position = behavior.timestamps(valid_speed_idx);
    x_pos = behavior.position.x(valid_speed_idx);
    y_pos = behavior.position.y(valid_speed_idx);
    fp_red = fp_red_interp(valid_speed_idx);
    fp_green = fp_green_interp(valid_speed_idx);
end

% Defining spatial bins
x_edges = linspace(0, 1, nBins+1);
y_edges = linspace(0, 1, nBins+1);


% Initialize maps

map_red = zeros(nBins,nBins);
map_green = zeros(nBins,nBins);
occupancy = zeros(nBins,nBins);
dt_sample = median(diff(timestamps_position));


% Fullfill maps
for i = 1:length(timestamps_position)
    x_idx = find(x_pos(i) >= x_edges(1:end-1) & x_pos(i) < x_edges(2:end));
    y_idx = find(y_pos(i) >= y_edges(1:end-1) & y_pos(i) < y_edges(2:end));

    if ~isempty(x_idx) && ~ isempty(y_idx)
        if isfield(fiber,'red')
            map_red(y_idx, x_idx) = map_red(y_idx, x_idx) + fp_red(i);
        end
        map_green(y_idx, x_idx) = map_green(y_idx, x_idx) + fp_green(i);
        occupancy(y_idx, x_idx) = occupancy(y_idx, x_idx) + dt_sample;
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
    red_normalization = imfilter(valid_mask_red, G, 'same');
    red_z = red_smoothed ./ red_normalization;
    red_z(isinf(red_z)) = NaN;
end

valid_mask_green = ~isnan(map_green_avg);
map_green_avg_zero = map_green_avg;
map_green_avg_zero(~valid_mask_green) = 0;
green_smoothed = imfilter(map_green_avg_zero, G, 'same');
green_normalization = imfilter(valid_mask_green, G, 'same');
green_z = green_smoothed ./ green_normalization;
green_z(isinf(green_z)) = NaN;



%% restructure into cell info data type

session = loadSession;

fiberMaps = [];

fiberMaps.params.smooth = sigma;
fiberMaps.params.minTime = minTime;
fiberMaps.params.nBins = nBins;
fiberMaps.params.speedThresh = speedThresh;

fiberMaps.green.occupancy = occupancy;
fiberMaps.green.count = map_green;
fiberMaps.green.z = green_z;

if isfield(fiber,'red')
    fiberMaps.red.occupancy = occupancy;
    fiberMaps.red.count = map_red;
    fiberMaps.red.z = red_z;
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