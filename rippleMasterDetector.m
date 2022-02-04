function [ripples] = rippleMasterDetector(varargin)
%   rippleMasterDetector - Wrapped function to compute different
%                           characteristics about hippocampal ripples (100 ~ 200 Hz oscillations).
%
% USAGE
%   [ripples] = rippleMasterDetector(lfp.data,lfp.timestamps,<options>)
%       OR
%   [ripples] = rippleMasterDetector(basepath,channel,<options>) Probably
%      
%   Probably needs to be changed the usage
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria
%
% INPUTS - note these are NOT name-value pairs... just raw values
%    lfp            unfiltered LFP (one channel) to use
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see tables below)
%
%    OR
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection (0-indexed, a la neuroscope)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100]). 
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization (default = all)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection 
%                   (default = [130 200])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'plotType'   1=original version (several plots); 2=only raw lfp
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps 
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector
% SEE ALSO
%
%    See also bz_Filter, bz_RippleStats, bz_SaveRippleEvents, bz_PlotRippleStats.
%   
%   Develop by Manu Valero and Pablo Abad 2022. Buzsaki Lab.
warning('this function is under development and may not work... yet')

% Default values
p = inputParser;
addParameter(p,'thresholds',[2 5],@isnumeric)
addParameter(p,'durations',[30 100],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'frequency',1250,@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'show','off',@isstr)
addParameter(p,'noise',[],@ismatrix)
addParameter(p,'passband',[130 200],@isnumeric)
addParameter(p,'EMGThresh',.9,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'minDuration',20,@isnumeric);
addParameter(p,'plotType',2,@isnumeric);

if isstr(varargin{1})  % if first arg is basepath
    addRequired(p,'basepath',@isstr)
    addRequired(p,'channel',@isnumeric)    
    parse(p,varargin{:})
    basename = bz_BasenameFromBasepath(p.Results.basepath);
    basepath = p.Results.basepath;
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    lfp = bz_GetLFP(p.Results.channel-1,'basepath',p.Results.basepath,'basename',basename);%currently cannot take path inputs
    signal = bz_Filter(double(lfp.data),'filter','butter','passband',passband,'order', 3);
    timestamps = lfp.timestamps;
    channel = p.Results.channel;
elseif isnumeric(varargin{1}) % if first arg is filtered LFP
    addRequired(p,'lfp',@isnumeric)
    addRequired(p,'timestamps',@isnumeric)
    parse(p,varargin{:})
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    signal = bz_Filter(double(p.Results.lfp),'filter','butter','passband',passband,'order', 3);
    timestamps = p.Results.timestamps;
    basepath = pwd;
    basename = bz_BasenameFromBasepath(basepath);
elseif isstruct(varargin{1})
    % Added by Pablo Abad to manage when first input if lfp file, not
    % filtered
    addRequired(p,'lfp',@isstruct);
    parse(p,varargin{:})
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    timestamps = p.Results.lfp.timestamps;
    EMGThres  = p.Results.EMGThresh;
    signal = bz_Filter(double(p.Results.lfp.data),'filter','butter','passband',passband,'order',3);
    basepath = pwd;
    basename = bz_BasenameFromBasepath(basepath);
    
end

frequency = p.Results.frequency;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
thresholds = p.Results.thresholds;
durations = p.Results.durations;
minRippleDuration = p.Results.minDuration;
plotType = p.Results.plotType;


%% Load session metadata
session = sessionTemplate(basepath,'showGUI',false);

% %% Load best Ripple Channel based on hippocampalLayers
try
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        disp('Hippocampal layers file found. Loading file !');
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
        rippleChannel = hippocampalLayers.layers{hippocampalLayers.bestShank}.pyramidal;
    end
catch
    rippleChannel = [];
end
% Just plotting purposes
plotRippleChannel_temp('rippleChannel',channel);
ripples = bz_FindRipples(basepath,rippleChannel,'thresholds',thresholds,'passband',passband,...
    'EMGThresh',EMGThresh,'durations',durations, 'saveMat',true);
ripples = removeArtifactsFromEvents(ripples);
ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',2);

%% Stats
lfp = bz_GetLFP('all');
filtered = bz_Filter(lfp,'channels',channel-1,'filter','butter','passband',passband,'order',3);
[maps,data,stats] = bz_RippleStats(filtered.data,filtered.timestamps,ripples);
ripples.maps = maps;
ripples.data = data;
ripples.stats = stats;
bz_PlotRippleStats(ripples.maps, ripples.data, ripples.stats);


if saveMat
    disp('Saving Results...');
    save([session.general.name , '.ripples.events.mat']);
end


end

