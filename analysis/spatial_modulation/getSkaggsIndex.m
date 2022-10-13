function [skaggs] = getSkaggsIndex(varargin)
%
%   Computes skaggs index.
%
%   USAGE
%      skaggs = getSkaggsIndex(<options>); 
%
%   INPUT
%       z :              firingMap
%       z_UnSmooth :     firingMap UnSmooth
%       occupancy :      position matrix
%       minTime :        minimum time in a bin to be included for analysis
%       unsmoothed_analysis : include or not unsmooth firing map. Default:
%                               false
%   OUTPUT
%
%       skaggs 
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'z_UnSmooth',[]);
addParameter(p,'occupancyUnSmooth',[]);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'unsmoothed_analysis',false,@islogical);
addParameter(p,'time',[]);

parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
z_UnSmooth = p.Results.z_UnSmooth;
occupancy = p.Results.occupancyUnSmooth;
minTime = p.Results.minTime;
unsmoothed_analysis = p.Results.unsmoothed_analysis;
time = p.Results.time;

skaggs = [];
skaggs.bitsPerSec = NaN;
skaggs.bitsPerSpike = NaN;

% We need to get the real duration of the recording and the real number of
% spikes



%% Compute Skaggs
nanmask = (z > 0) & time >= minTime;


duration = sum(occupancy(nanmask));
meanRate = sum(z(nanmask).*occupancy(nanmask))/duration;
if meanRate > 0.0
    p_x = occupancy ./ duration;
    p_r = z ./ meanRate;
    dummy = p_x .* z;
    idx = dummy > 0;
    skaggs.bitsPerSec = sum(dummy(idx).*log2(p_r(idx)));
    skaggs.bitsPerSpike = skaggs.bitsPerSec / meanRate;
end

    
end
