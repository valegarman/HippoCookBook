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
%   SharpWaves are detected based on the detected ripples. Radiatum lfp
%   signal is filtered (default [2 10] Hz) and szcore of the signal is
%   computed. SW peak is detected as the time when zscore rad signal
%   exceeds SWthreshold(2) during the ocurrence of a ripple (SW.peaks). If a ripple
%   does not have an associated SharpWave, nan values are in play. Onset
%   and offset of the sharpwave is computed as the time when radiatum
%   zscore signal first crosses threshold(1) both before and after the
%   peaks ( SW.timestamps)
%
%
%
%
%
%
%
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
addParameter(p,'SWthresholds',[-0.5 -2], @isnumeric);
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
addParameter(p,'rippleStats',false,@islogical);
addParameter(p,'debug',false,@islogical);
addParameter(p,'eventSpikeThreshold',1,@isnumeric);

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

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing Ripples
%%%%%%%%%%%%%%%%%%%%%%%%
ripples = findRipples(rippleChannel,'thresholds',thresholds,'passband',passband,...
    'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false);
ripples = removeArtifactsFromEvents(ripples);
ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',eventSpikeThreshold);

%% Ripple Stats
if rippleStats
    lfp = bz_GetLFP('all');
    filtered = bz_Filter(lfp,'channels',rippleChannel,'filter','butter','passband',passband,'order',3);
    [maps,data,stats] = bz_RippleStats(filtered.data,filtered.timestamps,ripples);
    ripples.maps = maps;
    ripples.data = data;
    ripples.stats = stats;
    bz_PlotRippleStats(ripples.maps, ripples.data, ripples.stats);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing SharpWaves
%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize of SW
SW = [];
SW.detectorinfo = [];
SW.timestamps = nan(size(ripples.timestamps,1),size(ripples.timestamps,2));
SW.peaks = nan(size(ripples.peaks,1),1);
SW.peakZScore = nan(size(ripples.peaks,1),1);

% Getting zscore signal
filteredRipple = bz_Filter(getLFP(rippleChannel),'filter','butter','passband',passband,'order',3);
filteredSW = bz_Filter(getLFP(SWChannel),'filter','butter','passband',SWpassband,'order',3);
zRipple = zscore(filteredRipple.data);
zSW = zscore(filteredSW.data);
ts = filteredRipple.timestamps;

eliminatedSW = false;
for i = 1:size(ripples.timestamps,1)
    dur = 0.01;
    signalSW = zSW(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
    signalR = zRipple(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
    t1SW = ts(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
    negPeak = find(signalSW < SWthresholds(2));

    if ~isempty(negPeak)
        [minP,indP] = min(signalSW);
        timePeak1 = t1SW(indP);
        SW.peaks(i) = t1SW(indP);
        SW.peakZScore(i) = signalSW(indP);
        disp('SharpWaveDetected')
        dur = 0.1;
        posSignalSW = zSW(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
        tSW = ts(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp)); 
        
        [minPeak,indPeak] = min(posSignalSW);
        timePeak2 = tSW(indPeak);
        while timePeak1 ~= timePeak2
            posSignalSW(indPeak) = nan;
            [minPeak,indPeak] = min(posSignalSW);
            timePeak2 = tSW(indPeak);
        end
        % First value crossing threshold1 (first timestamp of SharpWave)
        % It can happens that no crosses the threshold
        if (any(posSignalSW(indPeak-1:-1:1) > SWthresholds(1)) && any(posSignalSW(indPeak+1:1:length(posSignalSW)) > SWthresholds(1)))
            eliminatedSW = false;
            for j = indPeak-1:-1:1
                if (posSignalSW(j) > SWthresholds(1))
                    posPeak1 = j;
                    SW.timestamps(i,1) = tSW(j);
                    break 
                end
            end
            % Second value crossing threshold1 (second timestamp of SharpWave)
            for j = indPeak+1:1:length(posSignalSW)
                if (posSignalSW(j) > SWthresholds(1))
                    posPeak2 = j;
                    SW.timestamps(i,2) = tSW(j);
                    break     
                end
            end  
        else
            % Any of before or after SW peak is not crossing the
            % threshold(1). We discard this SW
            eliminatedSW = true;
            SW.peaks(i) = nan;
        end
        if debug
            posSignalR = zRipple(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
            figure,
            subplot(2,1,1)
            plot(tSW,posSignalR,'k')
            xlim([tSW(1) tSW(end)])
            subplot(2,1,2)
            plot(tSW,posSignalSW,'k')
            hold on
            xlim([tSW(1) tSW(end)])
            ylim([-4 4])
            if ~eliminatedSW
                text(tSW(posPeak1),posSignalSW(posPeak1),'X','Color','g');
                text(tSW(posPeak2),posSignalSW(posPeak2),'X','Color','g');
                text(tSW(indPeak),posSignalSW(indPeak),'O','Color','g');
                title('SharpWave !!')
            else
                title('SharpWave but eliminated ..');
                text(tSW(posPeak1),posSignalSW(posPeak1),'X','Color','r');
                text(tSW(posPeak2),posSignalSW(posPeak2),'X','Color','r');
                text(tSW(indPeak),posSignalSW(indPeak),'O','Color','r');
            end
            pause
            close all
        end
        
    else
        disp('No SW')
        dur = 0.1;
        posSignalSW = zSW(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
        tSW = ts(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp)); 
        if debug
            posSignalR = zRipple(round((ripples.timestamps(i,1)-dur)*srLfp):round((ripples.timestamps(i,2)+dur)*srLfp));
            figure,
            subplot(2,1,1)
            plot(tSW,posSignalR,'k')
            xlim([tSW(1) tSW(end)])
            subplot(2,1,2)
            plot(tSW,posSignalSW,'k')
            hold on
            xlim([tSW(1) tSW(end)])
            ylim([-4 4])
            title('No SharpWave detected')
            pause
            close all
        end
    end
end


% % of Sharp Waves ocurring in the detected ripples
SWnumber = sum(~isnan(SW.timestamps(:,1)));
SWperc = (SWnumber / size(ripples.timestamps,1));

% Distance of peaks between ripples and SW
rippleSWdifference = ripples.peaks - SW.peaks;
SWbefore = find(rippleSWdifference > 0);
SWbeforePerc = length(SWbefore) / SWnumber;
SWafter = find(rippleSWdifference < 0);
SWafterPerc = length(SWafter) / SWnumber;
SWequal = find(rippleSWdifference == 0);
SWequalPerc = length(SWequal) / SWnumber;

%% OUTPUT SHARPWAVES
SW.stats.SWperc = SWperc;
SW.stats.SWbeforeperc = SWbeforePerc;
SW.stats.SWafterperc = SWafterPerc;
SW.stats.SWequalperc = SWequalPerc;

SW.detectorinfo.detectionparams = p.Results;
SW.detectorinfo.detectorname = 'rippleMasterDetector';
SW.detectorinfo.detectiondate = today;
SW.detectorinfo.detectionintervals = [];
SW.detectorinfo.detectionchannel = SWChannel;


if saveMat
    disp('Saving Ripples Results...');
    save([session.general.name , '.ripples.events.mat'],'ripples');
    
    disp('Saving SharpWaves Results...');
    save([session.general.name , '.sharpwaves.events.mat'],'SW');
end


end

