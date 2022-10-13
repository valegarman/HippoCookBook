function [map] = MapTint(positions,times,varargin)

% USAGE
% [firingMaps] = MapTint(positions,spikes,varargin)
% Calculates averaged firing map for a set of 2D positions based on tint
% algorithm
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - 
%   <options>      optional list of property-value pairs (see table below)
% ===================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'var2binby'       'position','direction','speed',or 'pxd' - the
%                       variables to be binned by, in order (the values at
%                       the centre of each bin to be out in grid_values):
%                       units of cm for position, degrees for direction,
%                       cm/s for speed - depend on pixels_per_metre in .pos
%                       file header * with y increasing upwards*.
%     'binSize'			[8], [8], [8 8] - the sizef of the bins in each
%                       grid, in order matching var2binby. Binsize units for 'position' are
%                       camera pixels, degrees for 'direction'. Square bins are assumed for
%                       'position', and the range of position is the area tracked by the
%                       camera (using window_min_x etc in .pos header)
%                       For 'direction' binsize should be a factor of 360.
%     'pos2use'			list of samples to be binned (e.g. the position
%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%                .cmBin                 cm/bins ratio
%
% Pablo Abad 2022. Based on bin_pos_data by DACQ2matlab_distrib

%% parse inputs
p=inputParser;
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'plt',true,@islogical);
addParameter(p,'binsize',[],@isnumeric);
addParameter(p,'bndbox',[],@isstruct);
addParameter(p,'smooth',7,@isnumeric);
addParameter(p,'var2binby','position',@ischar);
addParameter(p,'sample_rate',30,@isnumeric);
addParameter(p,'pixelsmetre',[],@isnumeric);
addParameter(p,'minTime',0,@isnumeric);


parse(p,varargin{:});

saveMat = p.Results.saveMat;
plt = p.Results.plt;
binsize = p.Results.binsize;
bndbox = p.Results.bndbox;
smooth = p.Results.smooth;
var2binby = p.Results.var2binby;
sample_rate = p.Results.sample_rate;
pixelsmetre = p.Results.pixelsmetre;
minTime = p.Results.minTime;

% Check range of data - relative units from top left are used in .pos (pos 'position' and 'pxd')
% values are in cm, need to be converted to pixels
win_max_x = (round(bndbox.xmax) * pixelsmetre) /100;
win_min_x = (round(bndbox.xmin) * pixelsmetre) /100;
win_max_y = (round(bndbox.ymax) * pixelsmetre) /100;
win_min_y = (round(bndbox.ymin) * pixelsmetre) /100;


% Some info about x , y and z
if isempty(positions) || size(positions,1) < 2
    return;
end

pointProcess = (isempty(positions) | size(positions,2) == 1);
t = positions(:,1);
x = positions(:,2);
if size(positions,2) >= 3
    y = positions(:,3);
else
    y = [];
end

% bin the position data
[pos_binned_array] = get_bin_pos_data(var2binby,binsize,(positions(:,2:3)*pixelsmetre)/100,1:size(positions,1),'bndbox',bndbox,'pixelsmetre',pixelsmetre);

% bin the Spike data
[clust_binned_array] = get_bin_pos_data(var2binby,binsize,(positions(:,2:3)*pixelsmetre)/100,times,'bndbox',bndbox,'pixelsmetre',pixelsmetre);

unsmoothed_rate = (clust_binned_array./pos_binned_array)*sample_rate;

% smooth the data
[smoothed_pos, smoothed_spikes, smoothed_rate] = smooth_field_plot(pos_binned_array,...
    clust_binned_array, smooth, 'boxcar');

smoothed_rate = smoothed_rate.*sample_rate;

peak_rate = max(max(smoothed_rate));
mean_rate = nanmean(nanmean(smoothed_rate));

% Computations
if isempty(y)
    % 1D
    map.x = win_min_x:binsize:win_max_x;
    map.count = smoothed_spikes;
    map.time = smoothed_pos;
    map.z = smoothed_rate;
    valid = map.time > minTime;
    map.countUnSmooth = clust_binned_array;
    map.timeUnSmooth = pos_binned_array;
    map.zUnSmooth = unsmoothed_rate;
else
    % 2D (x and y)
    map.x = win_min_x:binsize:win_max_x;
    map.y = win_min_y:binsize:win_max_y;
    map.count = smoothed_spikes;
    map.time = smoothed_pos;
    map.z = smoothed_rate;
    valid = map.time > minTime;
    map.countUnSmooth = clust_binned_array;
    map.timeUnSmooth = pos_binned_array; % to convert to seconds
    map.timeUnSmoothSec = map.timeUnSmooth/sample_rate;
    map.zUnSmooth = unsmoothed_rate;
end
try
    map.bndbox = bndbox;
    map.var2binby = var2binby;
    map.binsize = binsize;
    map.pixelsmetre = pixelsmetre;
catch
    
end
    
end
