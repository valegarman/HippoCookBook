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
%     'sr'   sampling rate (in Hz) (default = 1250Hz)
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
%     'detector'   'filter' (default), 'cnn', 'consensus'
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
%   Develop by Manu Valero and Pablo Abad 2022. Buzsaki Lab. Based on
%   bz_findRipples

warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'referenceChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thresholds',[1.5 3.5],@isnumeric);
addParameter(p,'durations',[30 100],@isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'sr',1250,@isnumeric);
addParameter(p,'passband',[120 200],@isnumeric);
addParameter(p,'hfo_passband',[200 500],@isnumeric);
addParameter(p,'EMGThresh',1,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minDuration',20,@isnumeric);
addParameter(p,'rippleStats',true,@islogical);
addParameter(p,'eventSpikeThreshold',.5);
addParameter(p,'eventSpikeThreshold_shanks','all');
addParameter(p,'force',false,@islogical);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'useCSD',false,@islogical);
addParameter(p,'stdThreshold',2,@isnumeric);
addParameter(p,'detector','filter',@ischar);
addParameter(p,'excludeIntervals',[],@isnumeric);
% -- cnn related --
addParameter(p,'cnn_channels', [], @isnumeric);
addParameter(p,'model_file', '', @ischar);
addParameter(p,'pred_every', 0.032, @isnumeric);
addParameter(p,'verbose', false, @islogical);
addParameter(p,'handle_overlap', 'mean', @ischar);
addParameter(p,'exec_env', '', @ischar);

parse(p,varargin{:})

basepath = p.Results.basepath;
rippleChannel = p.Results.rippleChannel;
referenceChannel = p.Results.referenceChannel;
SWChannel = p.Results.SWChannel;
thresholds = p.Results.thresholds;
durations = p.Results.durations;
restrict = p.Results.restrict;
sr = p.Results.sr;
passband = p.Results.passband;
hfo_passband = p.Results.hfo_passband;
EMGThresh = p.Results.EMGThresh;
saveMat = p.Results.saveMat;
minDuration = p.Results.minDuration;
rippleStats = p.Results.rippleStats;
eventSpikeThreshold = p.Results.eventSpikeThreshold;
eventSpikeThreshold_shanks = p.Results.eventSpikeThreshold_shanks;
force = p.Results.force;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
useCSD = p.Results.useCSD;
stdThreshold = p.Results.stdThreshold;
detector = p.Results.detector;
excludeIntervals = p.Results.excludeIntervals;
% -- cnn related --
cnn_channels = p.Results.cnn_channels;
pred_every = p.Results.pred_every;
verbose = p.Results.verbose;
model_file = p.Results.model_file;
handle_overlap = p.Results.handle_overlap;
exec_env = p.Results.exec_env;
dir_cnn = fileparts(which('detect_ripples_cnn.m'));


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

if skipStimulationPeriods && ~isempty(dir('*optogeneticPulses.events.mat'))
    targetFile = dir('*optogeneticPulses.events.mat'); load(targetFile.name);
    try
        restrict_temp = SubtractIntervals([0 Inf],optoPulses.stimulationEpochs);
        restrict =  ConsolidateIntervals([restrict; restrict_temp; restrict_temp]);
    catch
        restrict_temp = SubtractIntervals([0 Inf],pulses.stimulationEpochs);
        restrict =  ConsolidateIntervals([restrict; restrict_temp; restrict_temp]);
    end
end

if useCSD
    disp('Computing CSD...');
    rippleChannel = computeCSD([],'channels',rippleChannel);
end



%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detecting Ripples
%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(detector, 'filter')
    ripples = findRipples(rippleChannel,'thresholds',thresholds,'passband',passband,...
        'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false,'restrict',restrict,'frequency',sr,'excludeIntervals',excludeIntervals,'noise',referenceChannel, 'minDuration', minDuration);
elseif strcmp(detector, 'cnn')
    % Load channel configuration
    load([session.general.name,'.session.mat']);
    channelConf = session.extracellular.electrodeGroups.channels;
    % Get best shank
    if isvarname('hippocampalLayers')
        bestShank = hippocampalLayers.bestShank;
    elseif ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
        bestShank = hippocampalLayers.bestShank;
    else
        [hippocampalLayers] = getHippocampalLayers();
        bestShank = hippocampalLayers.bestShank;
    end
    % Get channels near pyr that fits in the shank
    channelsShank = channelConf{bestShank};
    rippleChPos = find(channelsShank==rippleChannel);
    if (rippleChPos < length(channelsShank)-5) && (rippleChPos > 2)
        cnn_channels = channelsShank(rippleChPos + [-2:5]);
    elseif (rippleChPos < length(channelsShank)-5)
        cnn_channels = channelsShank(1:8);
    elseif (rippleChPos > 2)
        cnn_channels = channelsShank(end-7:end);
    end
    % Get LFP
    lfp = getLFP(cnn_channels,'basepath',basepath);
    % Get model directory
    if isempty(model_file)
        model_file = fullfile(dir_cnn, 'cnn');
    end
    % Predict probability
    ripples = detect_ripples_cnn(double(lfp.data), lfp.samplingRate, 'model_file', model_file, ...
        'pred_every', pred_every, 'verbose', verbose, 'handle_overlap', handle_overlap, ...
        'exec_env', exec_env);
elseif strcmp(detector, 'filter2')
    % As in Castelli et al, 2025
    ripples = findRipples(rippleChannel,'thresholds',thresholds,'passband',passband,...
        'EMGThresh',0,'durations',[20 300], 'saveMat',false,'restrict',restrict,'frequency',sr,'excludeIntervals',excludeIntervals,'noise',[], 'minDuration', 0);
    lfp_ripple_reference =      getLFP([rippleChannel referenceChannel],'noPrompts', true);        % 
    lfp_ripple_reference_filt = bz_Filter(lfp_ripple_reference, 'passband', passband, 'filter', 'butter', 'order', 4);
    
    hfo_ripple_reference_filt = bz_Filter(lfp_ripple_reference, 'passband', hfo_passband, 'filter', 'butter', 'order', 4);

    % Use existing ripple events to calculate number of cycles and mean frequency
    numEvents = length(ripples.peaks);
    num_cycles = zeros(numEvents,1);        % num_cycles: number of cycles for each ripple
    mean_freq = zeros(numEvents,1);         % mean_freq: avarage frequency of each ripple

    phase = lfp_ripple_reference_filt.phase(:,1);
    
    for ii = 1:numEvents
        % Convert times to samples
        onset_sample = floor(ripples.timestamps(ii,1) * sr);
        offset_sample = ceil(ripples.timestamps(ii,2) * sr);
    
        % Get phase for event window
        event_phase = phase(onset_sample:offset_sample);
        
        % Unwrap phase to avoid discontinuities
        unwrapped_phase = unwrap(event_phase);
        
        % Number of cycles = total phase change / 2*pi
        phase_diff = unwrapped_phase(end) - unwrapped_phase(1);
        num_cycles(ii) = phase_diff/(2*pi);
        
        % Duration in seconds
        duration_s = (offset_sample - onset_sample)/sr;
        
        % Mean frequency of ripple event
        mean_freq(ii) = num_cycles(ii) / duration_s;

        % ripple power in ripple channel
        rip_seg = lfp_ripple_reference_filt.data(onset_sample:offset_sample,1);
        pow_rip_detect(ii) = mean(rip_seg.^2);

        % ripple power in reference channel
        rip_ref_seg = lfp_ripple_reference_filt.data(onset_sample:offset_sample,2);
        pow_rip_ref(ii) = mean(rip_ref_seg.^2);

        % HF power in reference channel
        hf_seg = hfo_ripple_reference_filt.data(onset_sample:offset_sample,1);
        pow_hf(ii) = mean(hf_seg.^2);
    end

    % 1. The ripple band power (derived from squaring the mean ripple amplitude) in the detection channel should exceed twice the magnitude obtained for the reference channel.  
    % 2. The mean frequency of the event should surpass 100 Hz; 
    % 3. The event must comprise a minimum of 4 complete ripple cycles; 
    % 4. The power in the ripple band should be at least double compared to the control high frequency band.

    criterion1 = pow_rip_detect > 2 * pow_rip_ref;
    criterion2 = mean_freq > 100;
    criterion3 = num_cycles >= 4;
    criterion4 = pow_rip_detect > 2 * pow_hf;
    
    validEvents = criterion1 & criterion2' & criterion3' & criterion4;
    
    valid_timestamps = ripples.timestamps(validEvents, :);
    valid_peaks = ripples.peaks(validEvents);
    
    % saving results
    ripples.timestamps = valid_timestamps;
    ripples.peaks = valid_peaks;
    ripples.peakNormedPower = ripples.peakNormedPower(validEvents);
    
    fprintf('valid events found: %d su %d\n', sum(validEvents), numEvents);

elseif strcmp(detector, 'consensus')
    
else
    warning([detector ' is not recognised as a detector, please chose between filter/cnn/consensus'])
    return
end

if skipStimulationPeriods
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

ripples = removeArtifactsFromEvents(ripples,'stdThreshold',stdThreshold);

if isnumeric(eventSpikeThreshold) || eventSpikeThreshold
    if islogical(eventSpikeThreshold)
        eventSpikeThreshold = 1;
    end
    ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',eventSpikeThreshold,'shanksID',eventSpikeThreshold_shanks);
end

% plotRippleChannel('rippleChannel',rippleChannel,'ripples',ripples); % to do, run this after ripple detection
% EventExplorer(pwd, ripples)

%% Ripple Stats
if rippleStats
    ripples = computeRippleStats('ripples',ripples,'rippleChannel',rippleChannel);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Computing SharpWaves
% %%%%%%%%%%%%%%%%%%%%%%%%%
% if
% try
%     SW = findSharpWaves('ripples',ripples,'rippleChannel',rippleChannel,'SWChannel',SWChannel,...
%         'passband',passband,'SWpassband',SWpassband);
% catch
%     disp('Not possible to compute Sharp Waves...');
% end
%% OUTPUT
if saveMat
    disp('Saving Ripples Results...');
    save([session.general.name , '.ripples.events.mat'],'ripples');
    % try
    %     disp('Saving SharpWaves Results...');
    %     save([session.general.name , '.sharpwaves.events.mat'],'SW');
    % catch
    % end
end


end
