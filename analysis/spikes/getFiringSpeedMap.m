function [speedMaps] = getFiringSpeedMap(positions,spikes,varargin)
% USAGE
% [speedMaps] = getFiringSpeedMap(positions,spikes,varargin)
% Calculates averaged firing map for a set of velocities 
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
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   speedMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%                .cmBin                 cm/bins ratio
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',0,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;


if isstruct(positions)
    x = positions.position.x;
    y = positions.position.y;
    ts = positions.timestamps;
end
  
%% Calculate
% Compute speed
[~, ~,~,vx, vy, ax, ay] = KalmanVel(x,y,ts,order);
speed = sqrt(vx.^2 + vy.^2);

positions = [ts speed];
% bin size: 2cm/s

% get firign rate maps
for unit = 1:length(spikes.times)
    map{unit} = Map(positions,spikes.times{unit},'nBins',nBins,'maxGap',maxGap,'minTime',minTime,'mode',mode,'maxDistance',maxDistance );
end

% cmBin = (max(positions{1}(:,2))-min(positions{1}(:,2)))/nBins;
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
speedMaps.UID = spikes.UID;
try speed.sessionName = spikes.sessionName;
catch
    speedMaps.sessionName = spikes.basename;
end
try
speedMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

speedMaps.params.smooth = smooth;
speedMaps.params.minTime = minTime;
speedMaps.params.nBins = nBins;
speedMaps.params.maxGap = maxGap;
speedMaps.params.mode = mode;
speedMaps.params.maxDistance = maxDistance;
% speedMaps.cmBin = cmBin;

for unit = 1:length(spikes.times)
    speedMaps.rateMaps{unit,1} = map{unit}.z;
    speedMaps.countMaps{unit,1} = map{unit}.count;
    speedMaps.occupancy{unit,1} = map{unit}.time;
end

if saveMat
   save([speedMaps.sessionName '.speedMapsAvg.cellinfo.mat'],'speedMaps'); 
end

end
