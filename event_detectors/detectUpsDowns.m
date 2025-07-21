function [UDStates] = detectUpsDowns(varargin)
% Detect up and down states from LFP.
%
% USAGE
%
%   [UDStates] = detectUD(varargin)
%
% INPUTS
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%                      whateverPath/baseName/
%                      folder must include files:
%                      baseName.lfp
%                      baseName.sessionInfo.mat -or- baseName.xml
%                      recommended (will prompt if not included unless 'noPrompts')
%   ch             - Detection channel, by default choose the channel with the
%                      highest gamma/delta anticorrelation and highest gamma
%                      std.
%   noPrompts      - true/false disable any user prompts (default: true)
%   NREMInts       - Interval of times for NREM (seconds) 
%                      (Default: loaded from SleepState.states.mat, 
%                      it uses all recording if doesnt exist).
%   pulPeriods     - Interval of times with stimulation pulses.
%                      (Default: load from OptogeneticPulses.events.mat, it uses all
%                      recording if doesn't exist).
%   filterparams   - Structure as filterparams.gamma and .delta with filter
%                       parameters (default, filterparams.gamma = [20 200]
%                       and filterparams.delta = [.5 8]).
%   smoothGamma    - Envelope gamma moving average window, in seconds (default .1)
%   gamm_thr       - Gamma evelope threshold, in std (default 1)
%   down_thr       - Down state threshold, in -std of the gamma envelope (default 0)
%   minEvDistance  - Minimum distance between up and down peaks, seconds (default .5)
%   deltaWaveThreshold 
%                  - Delta wave threshold of the lfp signal (default 0)
%   plotOpt        - true or false option for plots with Up and Down sizes 
%                       and averages (default false)
%   forceDetect    - true or false to force detection (avoid load previous 
%                       detection, default false)
%   saveMat        - true or false save a buzcode event file with results
%                       (default true)
%   spikeThreshold - scalar, if [] false, default .25
%   skipCluster    - vector with cell ID for non-hippocampal and down state active cells, default [] 
%
% OUTPUT
%   UDStates structure with the following fields:
%   ints.DOWN               - Down state intervals
%   ints.UP                 - UP state intervals
%   timestamps.DOWN         - Down state peaks.
%   timestamps.UP           - Up state peaks.
%   timestamps.deltaWave

%
%   MV-BuzsakiLab 2019
%
%
%   Updated by Pablo Abad 2022 to discard bad channels ( I think it was already done) from session
%   metadata and remove sessionInfo dependencies
%   Modified by Pablo Abad to take into account timestamps where lfp is
%   completely flat (2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Defaults and Parms
filterparams.gamma = [30 600]; % 20 400, 20 100 
filterparams.delta = [.5 8];

p = inputParser;
addParameter(p,'basepath',pwd,@isstring);
addParameter(p,'ch',[],@isnumeric);
addParameter(p,'NREMInts',[]);
addParameter(p,'pulPeriods',[]);
addParameter(p,'filterparams',filterparams,@isstruct);
addParameter(p,'smoothGamm',.1,@isnumeric);
addParameter(p,'gamm_thr',.5,@isnumeric);
addParameter(p,'down_thr',0,@isnumeric);
addParameter(p,'minEvDistance',.5,@isnumeric);
addParameter(p,'deltaWaveThreshold',0,@isnumeric);
addParameter(p,'plotOpt',false,@islogical);
addParameter(p,'forceDetect',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'noPrompts',true,@islogical);
addParameter(p,'spikeThreshold',.25,@isnumeric);
addParameter(p,'skipCluster',[],@isnumeric);
addParameter(p,'useparfor',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'skipStimulationPeriods',true,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
noPrompts = p.Results.noPrompts;
ch = p.Results.ch;
NREMInts = p.Results.NREMInts;
pulPeriods = p.Results.pulPeriods;
filterparams = p.Results.filterparams;
smoothGamm = p.Results.smoothGamm;
gamm_thr = p.Results.gamm_thr;
down_thr = p.Results.down_thr;
minEvDistance = p.Results.minEvDistance;
deltaWaveThreshold = p.Results.deltaWaveThreshold;
plotOpt = p.Results.plotOpt;
forceDetect = p.Results.forceDetect;
saveMat = p.Results.saveMat;
spikeThreshold = p.Results.spikeThreshold;
skipCluster = p.Results.skipCluster;
useparfor = p.Results.useparfor;
excludeIntervals = p.Results.excludeIntervals;
skipStimulationPeriods = p.Results.skipStimulationPeriods;

%% Collect pieces
if exist([basenameFromBasepath(pwd) '.UDStates.events.mat'],'file') && ~forceDetect
    disp('Up and down states already detected! Loading file.');
    load([basenameFromBasepath(pwd) '.UDStates.events.mat']);
    return
end
tlfp = getLFP(1,'basepath',basepath,'noPrompts',noPrompts);
if ~isempty(excludeIntervals)
    [a,b] = InIntervals(tlfp.timestamps,excludeIntervals);
%     tlfp.timestamps(a) = nan(1,length(find(a)));
    tlfp.data(a,:) = nan(1,length(find(a)));
end
if ischar(NREMInts) && strcmpi(NREMInts,'all')
    NREMInts = [0 tlfp.duration];                                          % Consider all file as NREMInts
elseif isempty(NREMInts)                                                   % if empty, try to generate
    [SleepState] = bz_LoadStates(basepath,'SleepState');    
    if ~isempty(SleepState)
        NREMInts = SleepState.ints.NREMstate;
    else
        NREMInts = [0 tlfp.duration];                                      % Consider all file as NREMInts
    end
elseif isnumeric(NREMInts)
    disp('Using NREM periods specified by user...');
else
    warning('NREM interval format not recognized!!');
end

if isempty(pulPeriods) && skipStimulationPeriods
    if ~isempty(dir('*.optogeneticPulses.events.mat'))
        f = dir('*.optogeneticPulses.events.mat');
        disp('Using stimulation periods from optogeneticPulses.events.mat file');
        load(f.name);
        try
            pulPeriods = optoPulses.stimulationEpochs;
        catch
            pulPeriods = pulses.stimulationEpochs;
        end
        
    elseif ~isempty(dir('.pulses.events.mat'))
        f = dir('*Pulses.events.mat');
        disp('Using stimulation periods from pulses.events.mat file');
        load(f.name);
        pulPeriods = pulses.intsPeriods;
    else
        warning('No pulses epochs detected!');
    end
end

% combine stimulation and NREM periods
indState = zeros(size(tlfp.timestamps));                                   % Exclude stimulations times for NREMInts
for ii = 1:size(NREMInts,1)                                                % NRMInts, asign ones to NREM
    indState(tlfp.timestamps >= NREMInts(ii,1) & tlfp.timestamps <= NREMInts(ii,2)) = 1;
end
for ii = 1:size(pulPeriods,1)                                              % NRMInts, asign zeros to stimulation times
    indState(tlfp.timestamps >= pulPeriods(ii,1) & tlfp.timestamps <= pulPeriods(ii,2)) = 0;
end                                                                    % start of pulses
NREMInts = [tlfp.timestamps(find(diff(indState)==1)+1) tlfp.timestamps(find(diff(indState)==-1))]; 
NREM_ts = find(indState); % NREM timestamps
clear tlfp

session = loadSession;
channels = 1:session.extracellular.nChannels;
if isempty(ch) % if no channel declared, pick channel with higher gamma std during NREM and out of stimulation periods %anticorrelation delta gamma
    disp('Selecting best channel... ');
    params.Fs = session.extracellular.srLfp; params.fpass = [1 400]; params.tapers = [3 5]; params.pad = 1;
    if useparfor
        parfor ii = 1:session.extracellular.nChannels
            if isfield(session.channelTags,'Bad') && any(channels(ii) == session.channelTags.Bad.channels)
                gdCorr(ii) = 0; stdGamma(ii) = 0;
            else
                fprintf(' **Channel %3.i/ %3.i, ',ii, session.extracellular.nChannels); %\n
                chlfp = getLFP(ii,'basepath',basepath,'noPrompts',noPrompts);
                [S,t,f] = mtspecgramc_fast(double(chlfp.data(NREM_ts)),[2 1],params);
                S(find(S == 0)) = NaN;
                S = 10 * log10(S)' + 60; 
                gamm = nansum(S(find(f >= filterparams.gamma(1) & f <= filterparams.gamma(2)),:));
                delt = nansum(S(find(f >= filterparams.delta(1) & f <= filterparams.delta(2)),:));
                avGamma(ii) = nanmean(gamm);
                stdGamma(ii) = nanstd(gamm);
                try gdCorr(ii) = corr(gamm',delt','Type','Spearman');
                catch gdCorr(ii) = 0;
                end
            end
        end 
    else
        for ii = 1:session.extracellular.nChannels
            if isfield(session.channelTags,'Bad') && any(channels(ii) == session.channelTags.Bad.channels)
                gdCorr(ii) = 0; stdGamma(ii) = 0;
            else
                fprintf(' **Channel %3.i/ %3.i, ',ii, session.extracellular.nChannels); %\n
                chlfp = getLFP(ii,'basepath',basepath,'noPrompts',noPrompts);
                [S,t,f] = mtspecgramc_fast(double(chlfp.data(NREM_ts)),[2 1],params);
                S(find(S == 0)) = NaN;
                S = 10 * log10(S)' + 60; 
                gamm = nansum(S(find(f >= filterparams.gamma(1) & f <= filterparams.gamma(2)),:));
                delt = nansum(S(find(f >= filterparams.delta(1) & f <= filterparams.delta(2)),:));
                avGamma(ii) = nanmean(gamm);
                stdGamma(ii) = nanstd(gamm);
                try gdCorr(ii) = corr(gamm',delt','Type','Spearman');
                catch gdCorr(ii) = 0;
                end
            end
        end
    end
    gdCorr(gdCorr==0) = 1;
    sdScore = stdGamma./gdCorr;
    sdScore(~isnan(sdScore)) = zscore(sdScore(~isnan(sdScore)));                                     % down score increase with std gamma and gamma delta anticorr
    sdScore(sdScore>4) = 0; % remove outlier
    [~,ch] = max(sdScore(sdScore<4));
    ch = channels(ch);
    fprintf(' Best channel: %3.i!! \n',ch);
end
clear gamm delt avGamma stdGamma gdCorr dscore

%% find states
lfp = getLFP(ch,'basepath',basepath,'noPrompts',noPrompts);
if ~isempty(excludeIntervals)
    [a,b] = InIntervals(lfp.timestamps,excludeIntervals);
%     lfp.timestamps(a) = [];
    lfp.data(a,:) = nan(1,length(find(a)));
end
gammaLFP = bz_Filter(lfp,'passband',filterparams.gamma,'filter','fir1','order',4);
smoothGamm = smoothGamm * session.extracellular.srLfp;
envGamm = movmean(gammaLFP.data.^2,smoothGamm);
% envGamm = bz_Filter(envGamm,'passband',filterparams.delta,'filter','fir1','order',1);
envGamm = (envGamm - mean(envGamm((NREM_ts))))/std(envGamm((NREM_ts)));    % normalize by NREM periods
% find UP STATES!!
[vPeaks,upPeaks] = findpeaks(envGamm,session.extracellular.srLfp,'MinPeakDistance',minEvDistance,'MinPeakHeight',gamm_thr);
% find DOWN STATES!!
[dPeaks,downValley] = findpeaks(-envGamm,session.extracellular.srLfp,'MinPeakDistance',minEvDistance,'MinPeakHeight',down_thr);
if ~isempty(deltaWaveThreshold)                                            % DELTA WAVE flag, down state above deltaWaveThreshold
    deltaLFP = bz_Filter(lfp,'passband',filterparams.delta,'filter','fir1','order',1);
    deltaLFPNorm = (deltaLFP.data - mean(deltaLFP.data(NREM_ts)))/std(deltaLFP.data(NREM_ts));
    downValley = downValley(find(deltaLFPNorm(int32(downValley * session.extracellular.srLfp)) > deltaWaveThreshold)); 
end
indDown = [];
for ii = 1:length(downValley)                                              % pair each DOWN state with a previous or subsequent UP state
    indDown(ii) = any(abs(downValley(ii) - upPeaks) < minEvDistance * 2);
end
downValley = downValley(find(indDown));
% discard downstates out of NREM_st periods
downValley = downValley(find(indState(int32(downValley * session.extracellular.srLfp))));
downValley(downValley<1 | downValley> lfp.timestamps(end)-1) = [];         % discard downstate if in the first or last second of recording
indUp = [];
for ii = 1:length(upPeaks)                                                 % pair each UP state with a previous or subsequent DOWN state
    indUp(ii) = any(abs(upPeaks(ii) - downValley) < minEvDistance * 2);
end
upPeaks = upPeaks(find(indUp));

%% get intervals and peaks
disp('Find intervals... ');
win = int32(minEvDistance/2 * session.extracellular.srLfp);
if ~exist('deltaLFPNorm')
    deltaLFP = bz_Filter(lfp,'passband',filterparams.delta,'filter','fir1','order',1);
    deltaLFPNorm = (deltaLFP.data - mean(deltaLFP.data(NREM_ts)))/std(deltaLFP.data(NREM_ts));
end
for ii = 1:length(downValley)                                              % find delta wave peaks
    dw = deltaLFPNorm(downValley(ii)*session.extracellular.srLfp-win:downValley(ii)*session.extracellular.srLfp+win);
    [~,dP] = max(dw);
    if dP == 1 || dP == length(dw)
        deltaPeak(ii) = nan;
    else
        deltaPeak(ii) = double(dP + downValley(ii)*session.extracellular.srLfp - win - 1)/ session.extracellular.srLfp;
    end
end

envGammFilt = bz_Filter(envGamm,'passband',filterparams.delta,'filter','fir1','order',1);
win = int32(minEvDistance * 2 * session.extracellular.srLfp);
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                       % Returns Zero-Crossing Indices Of Argument Vector
for ii = 1:length(downValley)                                              % find DOWN interval
    dw = envGammFilt(downValley(ii)*session.extracellular.srLfp-win:downValley(ii)*session.extracellular.srLfp+win);
    zx = zci(dw); zx = [zx; win; 1250*2];
    DOWN(ii,1) = double(max(zx(zx<=win)) + downValley(ii)*session.extracellular.srLfp - win - 1)/ session.extracellular.srLfp;
    DOWN(ii,2) = double(min(zx(zx>win)) + downValley(ii)*session.extracellular.srLfp - win - 1)/ session.extracellular.srLfp;
end

upPeaks(find(upPeaks<= double(win)/session.extracellular.srLfp)) = [];
upPeaks(find(upPeaks>= lfp.duration-double(win)/session.extracellular.srLfp)) = [];
for ii = 1:length(upPeaks)                                                 % find UP interval
    uw = envGammFilt(upPeaks(ii)*session.extracellular.srLfp-win:upPeaks(ii)*session.extracellular.srLfp+win);
    zx = zci(uw);
    UP(ii,1) = double(max(zx(zx<win)) + upPeaks(ii)*session.extracellular.srLfp - win - 1)/ session.extracellular.srLfp;
    if isempty(min(zx(zx>win)))
        UP(ii,2) = double(max(zx(zx<win)) + upPeaks(ii)*session.extracellular.srLfp - win - 1)/ session.extracellular.srLfp + 0.1;
    end
end

%% Spike thresholding
if ~isempty(spikeThreshold)
    downWinSize = 0.4;
    spikes = loadSpikes;
    if ~isempty(skipCluster)                                               % skip non-cortical cluster and DownStateActive cells
        for ii = 1:length(skipCluster)
            spikes.times{skipCluster(ii)} = [];
        end
    end
    
    allspikes = sort(cat(1,spikes.times{:}));
    spikemat = bz_SpktToSpkmat(spikes.times,'binsize',0.03,'overlap',6);
    sSpkMat = sum(spikemat.data,2)/size(spikemat.data,2);
    downValley(downValley<downWinSize) = [];
    
    disp('Population cell response... ');
    downPopResponse = [];
    tic
    if useparfor
        parfor (ii = 1: length(downValley),6)
            try a = sSpkMat(spikemat.timestamps>=downValley(ii)-downWinSize ...
                & spikemat.timestamps<=downValley(ii)+downWinSize);
                downPopResponse(ii,:) = zscore(a(int32(1:downWinSize*2/(mean(diff(spikemat.timestamps)))-1)));
            end
        end
        toc
    else
         for ii = 1: length(downValley)
            try a = sSpkMat(spikemat.timestamps>=downValley(ii)-downWinSize ...
                & spikemat.timestamps<=downValley(ii)+downWinSize);
                downPopResponse(ii,:) = zscore(a(int32(1:downWinSize*2/(mean(diff(spikemat.timestamps)))-1)));
            end
         end
    end
    t_ds = linspace(-downWinSize,downWinSize,size(downPopResponse,2));
    
    % sort response according to population response
    [~,idx]=sort(mean(downPopResponse(:,t_ds>-downWinSize/8 & t_ds<downWinSize/8),2));
    
    if plotOpt
        figure
        imagesc(t_ds,1:length(downValley),downPopResponse(idx,:),[-3 3]);
        hold on
        plot([t_ds(1) t_ds(end)],[length(downValley)*spikeThreshold length(downValley)*spikeThreshold],'w','LineWidth',2)
        xlabel('s'); ylabel('# events');
        mkdir('SummaryFigures'); % create folder    
        saveas(gcf,'SummaryFigures\UDStates.png');
    end
    
    classEvs = sort(idx(1:int32(length(downValley)*spikeThreshold)));
end
%% Saving...

UDStates.ints.DOWN = DOWN(classEvs,:);
UDStates.ints.UP = UP;
UDStates.timestamps.DOWN = downValley(classEvs);
UDStates.timestamps.UP = upPeaks;
UDStates.timestamps.deltaWave = deltaPeak(~isnan(deltaPeak(classEvs)))';
UDStates.detectionInfo.ch = ch;
UDStates.detectionInfo.spikeThreshold = spikeThreshold;
UDStates.detectionInfo.gamm_thr = gamm_thr;
UDStates.detectionInfo.down_thr = down_thr;
UDStates.detectionInfo.deltaWaveThreshold = deltaWaveThreshold;
if saveMat
    save([session.general.name '.UDStates.events.mat'],'UDStates');
end

% %% check parameters
%t_in = [28*60+44.183 28*60+49.183];
% t_in = [1713 1718];
% 
% figure
% subplot(3,1,1)
% plot(lfp.timestamps,lfp.data)
% xlim([t_in]);
% subplot(3,1,2)
% hold on
% plot(lfp.timestamps,envGamm)
% plot(downValley,zeros(size(downValley)),'or')
% plot(upPeaks,zeros(size(upPeaks)),'og')
% xlim([t_in]);


% %% plots
% if plotOpt
%     % size histograms
%     xbins = [0:0.04:1.5];
%     figure
%     hold on
%     upH = histogram(diff(UP'),xbins,'FaceColor',[.3 .4 .8],'EdgeColor','none','FaceAlpha',.5);
%     doH = histogram(diff(DOWN'),xbins,'FaceColor',[.8 .4 .3],'EdgeColor','none','FaceAlpha',.5);
%     legend([upH doH], 'Up', 'Down');
%     ylabel('#'); xlabel('s');
%     
%     % all down
%     win = int32(minEvDistance * sessionInfo.rates.lfp);
%     for ii = 1:length(downValley)
%         AD(:,ii) = lfp.data(downValley(ii)*lfp.samplingRate-win : downValley(ii)*lfp.samplingRate+win);
%     end
%     xt = linspace(-minEvDistance,minEvDistance,size(AD,1));
%     [~,idxD] = sort(AD(win+1,:));
%     % all UP
%     for ii = 1:length(upPeaks)
%         AU(:,ii) = lfp.data(upPeaks(ii)*lfp.samplingRate-win : upPeaks(ii)*lfp.samplingRate+win);
%     end
%     [~,idxU] = sort(AU(win+1,:));
%     
%     figure
%     subplot(3,2,1)
%     plot(xt,mean(AD,2)/1000);
%     set(gca,'XTick',[]); ylabel('mV'); title('Down','FontWeight','normal');
%     subplot(3,2,2)
%     plot(xt,mean(AU,2)/1000);
%     set(gca,'XTick',[]); ylabel('mV'); title('Up','FontWeight','normal');
%     subplot(3,2,[3 5])
%     imagesc(xt,1:size(AD,2),AD(:,idxD)',[-abs(max(AD(:))) abs(max(AD(:)))]);
%     colormap jet
%     ylabel('#'); xlabel('s');
%     subplot(3,2,[4 6])
%     imagesc(xt,1:size(AU,2),AU(:,idxU)',[-abs(max(AU(:))) abs(max(AU(:)))]);
%     xlabel('s');
% end

end


