function [fiberMaps] = getFiberLinearSpatialMap(behavior,fiber,varargin)

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
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
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
    positions = behavior.maps;
    % maybe behavior.lin (to take all the maze and not diciding two
    % halves)
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

    timestamps_position{iCond} = positions{iCond}(:,1);
    x_pos{iCond} = positions{iCond}(:,2);

    % We need to separate between left and right trials
    idx = ismember(behavior.timestamps,timestamps_position{iCond});
    fp_red{iCond} = fp_red_interp(idx);
    fp_green{iCond} = fp_green_interp(idx);

    valid_speed_idx = v > velocity_thresh;
    
    timestamps_position{iCond} = timestamps_position{iCond}(valid_speed_idx);
    x_pos{iCond} = x_pos{iCond}(valid_speed_idx);        
    fp_red{iCond} = fp_red{iCond}(valid_speed_idx);
    fp_green{iCond} = fp_green{iCond}(valid_speed_idx);

end

% Defining spatial bins

edges = linspace(round(min(behavior.position.lin)),round(max(behavior.position.lin)), nBins+1);

% Fullfill maps
for j = 1:size(positions,2)
    % Initialize maps
    map_red{j} = zeros(1,nBins);
    map_green{j} = zeros(1,nBins);
    occupancy{j} = zeros(1,nBins);
    dt_sample{j} = median(diff(timestamps_position{j}));


    for i = 1:length(timestamps_position{j})
        x_idx = find(x_pos{j}(i) >= edges(1:end-1) & x_pos{j}(i) < edges(2:end));
    
        if ~isempty(x_idx)
            if isfield(fiber,'red')
                map_red{j}(x_idx) = map_red{j}(x_idx) + fp_red{j}(i);
            end
            map_green{j}( x_idx) = map_green{j}(x_idx) + fp_green{j}(i);
            occupancy{j}(x_idx) = occupancy{j}(x_idx) + dt_sample{j};

        end
    end

    % Compute mean signal per bin
    occupancy{j}(occupancy{j} < min_occupancy) = NaN;
    if normalize_by_occupancy
        map_red_avg{j} = map_red{j} ./ occupancy{j};
        map_green_avg{j} = map_green{j} ./ occupancy{j};
    else
        map_red_avg{j} = map_red{j};
        map_green_avg{j} = map_green{j};
    end

    % Smnoothing
    kernelSize = 7;  % impar (≈ 3*sigma)
    x = -(kernelSize-1)/2:(kernelSize-1)/2;
    G = exp(-(x.^2)/(2*sigma^2));
    G = G / sum(G);   % normalizar

    if isfield(fiber,'red')
        valid_mask_red{j} = ~isnan(map_red_avg{j});
        map_red_avg_zero{j} = map_red_avg{j};
        map_red_avg_zero{j}(~valid_mask_red{j}) = 0;

        red_smoothed{j} = imfilter(map_red_avg_zero{j}, G, 'same');
        red_normalization{j} = imfilter(valid_mask_red{j}, G, 'same');
        red_z{j} = red_smoothed {j} ./ red_normalization{j};
    end
    
    valid_mask_green{j} = ~isnan(map_green_avg{j});
    map_green_avg_zero{j} = map_green_avg{j};
    map_green_avg_zero{j}(~valid_mask_green{j}) = 0;
    green_smoothed{j} = imfilter(map_green_avg_zero{j}, G, 'same');
    green_normalization{j} = imfilter(valid_mask_green{j}, G, 'same');
    green_z{j} = green_smoothed {j} ./ green_normalization{j};
end

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
   save([session.general.name '.fiberLinearMapAvg.mat'],'fiberMaps'); 
end

%% Plotting

if plt

    xtrack = linspace(round(min(behavior.position.lin)),round(max(behavior.position.lin)),nBins);
    for j = 1:length(fiberMaps.green.z)
        % Green
        figure;
        plot(xtrack,fiberMaps.green.z{j},'k');
        title('Green fiber linear spatial map');
        xlabel('x(cm)');
        ylabel('dFF');
        ylim([0 1])
    
        if saveFig
            saveas(gca,['SummaryFigures\fiber_green_linear_spatial_map_',num2str(j),'.png'])
        end
        
        % Red
        if isfield(fiber,'red')
            figure;
            plot(xtrack,fiberMaps.red.z{j},'k');
            title('Red fiber linear spatial map');
            xlabel('X bin');
            ylabel('Y bin');
            ylim([0 1])
        end
    
        if saveFig
            saveas(gca,['SummaryFigures\fiber_red_linear_spatial_map_',num2str(j),'.png'])
        end
    end
end


end