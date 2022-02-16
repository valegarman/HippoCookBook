function [ripples] = rippleMasterDetector(varargin)
%   rippleMasterDetector - Wrapped function to compute different
%                           characteristics about hippocampal ripples (100 ~ 200 Hz oscillations).
%
% USAGE
%   [ripples] = rippleMasterDetector(<options>)
%   
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
%    <options>      optional list of property-value pairs (see tables below)
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

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thresholds',[2 5],@isnumeric);
addParameter(p,'SWthresholds',[0.5 -2], @isnumeric);
addParameter(p,'durations',[20 250],@isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'stdev',[],@isnumeric);
addParameter(p,'show','off',@isstr);
addParameter(p,'noise',[],@ismatrix);
addParameter(p,'passband',[80 200],@isnumeric);
addParameter(p,'SWpassband',[2 10],@isnumeric);
addParameter(p,'EMGThresh',1,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minDuration',20,@isnumeric);
addParameter(p,'plotType',2,@isnumeric);
addParameter(p,'srLfp',1250,@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thresholds = p.Results.thresholds;
SWthresholds = p.Results.SWthresholds;
durations = p.Results.durations;
restrict = p.Results.restrict;
frequency = p.Results.frequency;
stdev = p.Results.stdev;
show = p.Results.show;
noise = p.Results.noise;
passband = p.Results.passband;
SWpassband = p.Results.SWpassband;
EMGThresh = p.Results.EMGThresh;
saveMat = p.Results.saveMat;
minDuration = p.Results.minDuration;
plotType = p.Results.plotType;
srLfp = p.Results.srLfp;
%% Load Session Metadata and several variables if not provided
session = sessionTemplate(basepath,'showGUI',false);

% Ripple and SW Channel are loaded separately in case we want to provide
% only one of the
if isempty(rippleChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    rippleChannel = hippocampalLayers.layers{hippocampalLayers.bestShank}.pyramidal;
end

if isempty(SWChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    SWChannel = hippocampalLayers.layers{hippocampalLayers.bestShank}.radiatum;
end


% Just plotting purposes
plotRippleChannel('rippleChannel',rippleChannel);

% Computing Ripples
ripples = findRipples(rippleChannel,'thresholds',thresholds,'passband',passband,...
    'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false);
ripples = removeArtifactsFromEvents(ripples);
ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',2);

% Computing SharpWaves
lfpRipple = getLFP(rippleChannel);
lfpSW = getLFP(SWChannel);

filteredRipple = bz_Filter(double(lfpRipple.data),'filter','butter','passband',passband,'order',3);
filteredSW = bz_Filter(double(lfpSW.data),'filter','butter','passband',SWpassband,'order',3);

zRipple = zscore(filteredRipple);
zSW = zscore(filteredSW);

for i = 1:size(ripples.timestamps,1)
    signal = zSW(ripples.timestamps(i,1)*srLfp:ripples.timestamps(i,2)*srLfp);
    negPeak = find(zSW(signal < SWthresholds(2));
    if ~isempty(negPeak)
        disp('SharpWaveDetected')
        posPeak1 = find(signal < SWthresholds(1),1);
        posPeak2 = find(signal < SWthresholds(1),end);
    end
end




%% Stats
lfp = bz_GetLFP('all');
filtered = bz_Filter(lfp,'channels',channel,'filter','butter','passband',passband,'order',3);
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

