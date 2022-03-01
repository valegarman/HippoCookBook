function [SW] = findSharpWaves(varargin)
%   findSharpWaves - SharpWaves are detected based on the detected ripples.
%   Radiatum lfp signal is filtered (default [2 10] Hz) and szcore of the signal is
%   computed. SW peak is detected as the time when zscore rad signal
%   exceeds SWthreshold(2) during the ocurrence of a ripple (SW.peaks). If a ripple
%   does not have an associated SharpWave, nan values are in play. Onset
%   and offset of the sharpwave is computed as the time when radiatum
%   zscore signal first crosses threshold(1) both before and after the
%   peaks ( SW.timestamps)
%
% USAGE
%   [SW] = findSharpWaves(<options>)
%
% INPUTS
%   ripples - ripples struct
%   rippleChannel - channel for ripple detection
%   SWChannel - Channel for SW detection
%   passband - frequency range for ripple 
%   SWpassband - frequency range for sharp waves
%   debug - logical for debugging, (default false)
%   
%
% OUTPUT
% SW            buzcode format .event. struct with the followinf fields
%               .timestamps
%               .peaks
%               .detectorName
%               .stats
%
% Developed by Pablo Abad and Manuel Valero 2022. Buzsaki Lab.
%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ripples',[],@isstruct);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'passband',[120 200], @isnumeric);
addParameter(p,'SWpassband',[2 10], @isnumeric);
addParameter(p,'SWthresholds',[-0.5 -2],@isnumeric);
addParameter(p,'debug',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
ripples = p.Results.ripples;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
passband = p.Results.passband;
SWpassband = p.Results.SWpassband;
SWthresholds = p.Results.SWthresholds;
debug = p.Results.debug;
%% Load Session Metadata
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
srLfp = session.extracellular.srLfp;
%% Initialization of SW
SW = [];
SW.detectorinfo = [];
SW.timestamps = nan(size(ripples.timestamps,1),size(ripples.timestamps,2));
SW.peaks = nan(size(ripples.peaks,1),1);
SW.peakZScore = nan(size(ripples.peaks,1),1);

%% Computing Sharp Waves
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
%         disp('SharpWaveDetected')
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
%         disp('No SW')
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
SW.detectorinfo.detectiondate = date;
SW.detectorinfo.detectionintervals = [];
SW.detectorinfo.detectionchannel = SWChannel;

end

