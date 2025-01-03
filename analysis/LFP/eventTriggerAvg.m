function [psth] = eventTriggerAvg(timestamps,varargin)
% Computes event-triggered average for timestamps entered as inputs.
% USAGE
%   [psth] = eventTriggeredAvg(timestamps,<options>)
%
% INPUTS
%   timestamps - mx2 matrix indicating timetamps (in seconds) over which to
%                   compute psth
%   
% <OPTIONALS>
%   basepath - default pwd
%   lfp - buzcode spikes structure
%   numRep - For bootstraping, default 500. If 0, no bootstraping
%   binSize - In seconds, default 0.001 (1 ms)
%   winSize - In seconds, default 0.5 (500 ms)
%   rasterPlot - Default true
%   ratePlot - Default true
%   saveMat - default true
%   eventType - default, date, other options: ripples...
%   event_ints - interval around events timestamps to compute
%   baseline_ints - interval before events timestamps to compute baseline
%   minNumberOfPulses - minimum number of pulses to create pulses entry, default 100
%   win_Z - Interval arround events to compting Zscore mean and SD. By
%       default [-winSize/2 -events_ints(1)];
%
% OUTPUTS
%   eventAvg
%
% Developed by Pablo Abad 2024.
%% Defaults and Params
p = inputParser;

addRequired(p,'timestamps',@isnumeric);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'numRep',500,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',1,@isnumeric);
addParameter(p,'getRaster',true,@islogical);
addParameter(p,'ratePlot',true,@islogical);
addParameter(p,'winSizePlot',[-.1 .5],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'eventType',date,@ischar);
addParameter(p,'event_ints',[-0.02 0.02],@isnumeric);
addParameter(p,'baseline_ints',[-0.5 -0.46],@isnumeric);
addParameter(p,'minNumberOfPulses',100,@isnumeric);
% addParameter(p,'win_Z',[],@isnumeric);
addParameter(p,'restrictIntervals',[],@isnumeric);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'raster_time',[-0.250 0.250],@isnumeric);
addParameter(p,'salt_win',[0.01],@isscalar);
addParameter(p,'salt_binSize',[0.001],@isscalar);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','_psth',@ischar);
addParameter(p,'save_raster_as','_raster',@ischar);
addParameter(p,'sr',[],@isnumeric);


parse(p, timestamps,varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
getRaster = p.Results.getRaster;
ratePlot = p.Results.ratePlot;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
savePlot = p.Results.savePlot;
force = p.Results.force;
eventType = p.Results.eventType;
event_ints = p.Results.event_ints;
baseline_ints = p.Results.baseline_ints;
minNumberOfPulses = p.Results.minNumberOfPulses;
% win_Z = p.Results.win_Z;
restrictIntervals = p.Results.restrictIntervals;
bootsTrapCI = p.Results.bootsTrapCI;
salt_baseline = p.Results.salt_baseline;
raster_time = p.Results.raster_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;
save_raster_as = p.Results.save_raster_as;
sr = p.Results.sr;

%% Session Template
% Deal with inputs
prevPath = pwd;
cd(basepath);

if minNumberOfPulses < 2
    error('Number of pulses should be lager than 1');
end

session = loadSession;
if exist([session.general.name '.' eventType '_eventAvg.mat'],'file') ...
        && ~force
    disp(['Event average already computed for ', session.general.name, ' ', eventType,'. Loading file.']);
    load([session.general.name '.' eventType '_eventAvg.mat']);
    return
end

% if isempty(win_Z)
%     win_Z = [-winSize/2 -event_ints(1)];
% end

% default detection parameters
if strcmpi(eventType,'ripples')
    if isempty(timestamps)
        ripples = rippleMasterDetector;
        timestamps = ripples.peaks;
    end
    warning('Using default parameters for ripples!');
    binSize = 0.005;
    winSize = 1;
    winSizePlot = [-0.5 0.5];
    event_ints = [-0.025 0.025];
    baseline_ints = [-0.5 -0.5 + diff(event_ints)];
    % win_Z = [-0.5 -0.1];
end

%% LFP

if isempty(lfp)
    lfp = getLFP();

end


end
