function [ripples,SW] = rippleMasterDetector(varargin)
%   rippleMasterDetector - Wrapped function to compute different
%                           characteristics about hippocampal ripples (100
%                           ~ 200 Hz oscillations). It also computes
%                           SharpWaves based on the ripples detected.
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
%    exclusion criteria.
%
%   SharpWaves are detected based on the detected ripples by findSharpWaves. Radiatum lfp
%   signal is filtered (default [2 10] Hz) and szcore of the signal is
%   computed. SW peak is detected as the time when zscore rad signal
%   exceeds SWthreshold(2) during the ocurrence of a ripple (SW.peaks). If a ripple
%   does not have an associated SharpWave, nan values are in play. Onset
%   and offset of the sharpwave is computed as the time when radiatum
%   zscore signal first crosses threshold(1) both before and after the
%   peaks ( SW.timestamps)
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
%     'restrict'    interval used to compute normalization. By default,
%                       tries to remove stimulated periods. Otherwise, 'off'
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
%   SW              buzcode format .event. struct with the following fields
%                   .timestamps
%                   .detectorName
%                   .peaks
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
addParameter(p,'thresholds',[1.5 3.5],@isnumeric);
addParameter(p,'SWthresholds',[-0.5 -2], @isnumeric);
addParameter(p,'durations',[30 100],@isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'stdev',[],@isnumeric);
addParameter(p,'show','off',@isstr);
addParameter(p,'noise',[],@ismatrix);
addParameter(p,'passband',[120 200],@isnumeric);
addParameter(p,'SWpassband',[2 10],@isnumeric);
addParameter(p,'EMGThresh',1,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minDuration',20,@isnumeric);
addParameter(p,'plotType',2,@isnumeric);
addParameter(p,'srLfp',1250,@isnumeric);
addParameter(p,'rippleStats',true,@islogical);
addParameter(p,'debug',false,@islogical);
addParameter(p,'eventSpikeThreshold',1,@isnumeric);
addParameter(p,'force',false,@islogical);
addParameter(p,'removeRipplesStimulation',true,@islogical);
addParameter(p,'useCSD',false,@islogical);

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
rippleStats = p.Results.rippleStats;
debug = p.Results.debug;
eventSpikeThreshold = p.Results.eventSpikeThreshold;
force = p.Results.force;
removeRipplesStimulation = p.Results.removeRipplesStimulation;
useCSD = p.Results.useCSD;

%% Load Session Metadata and several variables if not provided
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);

if (exist([session.general.name '.ripples.events.mat'],'file') ...
        && ~force)
    disp(['Ripples already detected for ', session.general.name, '. Loading file.']);
    load([session.general.name '.ripples.events.mat']);
    return
end

% Ripple and SW Channel are loaded separately in case we want to provide
% only one of the
if isempty(rippleChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

if isempty(SWChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    SWChannel = hippocampalLayers.bestShankLayers.radiatum;
end

if removeRipplesStimulation
    targetFile = dir('*optogeneticPulses.events.mat'); load(targetFile.name);
    restrict_temp = SubtractIntervals([0 Inf],optoPulses.stimulationEpochs);
    restrict =  ConsolidateIntervals([restrict; restrict_temp]);
end

if useCSD
    disp('Computing CSD...');
    rippleChannel = computeCSD([],'channels',rippleChannel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing Ripples
%%%%%%%%%%%%%%%%%%%%%%%%
ripples = findRipples(rippleChannel,'thresholds',thresholds,'passband',passband,...
    'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false,'restrict',restrict,'frequency',frequency);
if removeRipplesStimulation
    try
        % Remove ripples durting stimulation artifacts
        if ~isempty(dir('*.optogeneticPulses.events.mat'))
            f = dir('*.optogeneticPulses.events.mat');
            disp('Using stimulation periods from optogeneticPulses.events.mat file');
            load(f.name);
            pulPeriods = optoPulses.stimulationEpochs;

        elseif ~isempty(dir('.pulses.events.mat'))
            f = dir('*Pulses.events.mat');
            disp('Using stimulation periods from pulses.events.mat file');
            load(f.name);
            pulPeriods = pulses.intsPeriods;
        else
            warning('No pulses epochs detected!');
        end
        nBefore = length(ripples.peaks);
        for i = 1:size(pulPeriods,1)
            a = InIntervals(ripples.peaks,pulPeriods(i,:));
            fieldsR = fields(ripples);
            for j = 1:size(fieldsR,1)
                if ~isstruct(ripples.(fieldsR{j})) && size(ripples.(fieldsR{j}),1) > 3
                    ripples.(fieldsR{j})(a,:) = [];
                end
            end
        end
        fprintf('%3.i/%3.i events during stimulation, discarted... \n',nBefore-length(ripples.peaks), nBefore); %\n
    catch
        warning('Not possible to remove ripples during stimulation epochs...');
    end
end
ripples = removeArtifactsFromEvents(ripples,'stdThreshold',1.5);
ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',eventSpikeThreshold);
plotRippleChannel('rippleChannel',rippleChannel,'ripples',ripples); % to do, run this after ripple detection
% EventExplorer(pwd, ripples)

%% Ripple Stats
if rippleStats
    ripples = computeRippleStats('ripples',ripples);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing SharpWaves
%%%%%%%%%%%%%%%%%%%%%%%%%
SW = findSharpWaves('ripples',ripples,'rippleChannel',rippleChannel,'SWChannel',SWChannel,...
    'passband',passband,'SWpassband',SWpassband);

%% OUTPUT
if saveMat
    disp('Saving Ripples Results...');
    save([session.general.name , '.ripples.events.mat'],'ripples');
    
    disp('Saving SharpWaves Results...');
    save([session.general.name , '.sharpwaves.events.mat'],'SW');
end


end
