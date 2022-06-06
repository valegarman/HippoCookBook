
function [thetaEpochs] = detectThetaEpochs(varargin)
% Detect theta/delta periods
% 
% INPUTS
% <optional>
% 'basepath'          Default pwd
% 'lfp'                 buzcode-formatted lfp structure (use bz_GetLFP)
%                           needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%                           If empty or no exist, look for lfp in basePath folder
% 'saveSummary'         Default true
% 'saveMat'             Detault true
% 'force'               Default false
% 'bandpass'            Default [6 12]
% 'powerThreshold'      Default 1 SD
% 'channel'             Numeric [ex, 5]; by default calls
%                           getHippocampalLayers and uses oriens.
% 'updateSleepStates'   Default true
% 'useCSD'              Default, true.
% 'discardRipples'      Discard ripples from nonTheta, default true.
% 
% OUTPUT
% thetaEpochs         states structure with theta epochs intervals
%
% Manu Valero 2022
% Include theta/delta ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'bandpass',[6 12], @isnumeric);
addParameter(p,'powerThreshold',[.8], @isnumeric);
addParameter(p,'channel',[],@isnumeric);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'updateSleepStates',true,@islogical);
addParameter(p,'useCSD',true,@islogical);
addParameter(p,'delta_bandpass',[1 3],@isnumeric);
addParameter(p,'use_ratioThetaDelta',true,@islogical);
addParameter(p,'threshold_noise',2.5,@isnumeric);
addParameter(p,'uselog10Power',true,@islogical);
addParameter(p,'powerThreshold_nonTheta',.5,@isnumeric);
addParameter(p,'discardRipples',true,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
theta_bandpass = p.Results.bandpass;
powerThresh = p.Results.powerThreshold;
channel = p.Results.channel;
plotting = p.Results.plotting;
updateSleepStates = p.Results.updateSleepStates;
useCSD = p.Results.useCSD;
delta_bandpass = p.Results.delta_bandpass;
use_ratioThetaDelta = p.Results.use_ratioThetaDelta;
threshold_noise =p.Results.threshold_noise;
uselog10Power = p.Results.uselog10Power;
powerThreshold_nonTheta = p.Results.powerThreshold_nonTheta;
discardRipples = p.Results.discardRipples;

% Deal with inputs
prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.thetaEpochs.states.mat');
if ~isempty(targetFile) && ~force
    disp('Theta epochs already detected! Loading file.');
    load(targetFile.name);
    return
end

% 
if isempty(channel)
    hippocampalLayers = getHippocampalLayers;
    channel = hippocampalLayers.bestShankLayers.slm;
end

if isempty(lfp) && ~useCSD
    lfpT = getLFP(channel,'noPrompts',true);
elseif useCSD
    disp('Computing CSD...');
    lfpT = computeCSD(lfp,'channels',channel);
else
    warning('CSD estimation not possible. Using LFP...');
end

samplingRate = lfpT.samplingRate;
[wave,f,t,~,wphases,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,theta_bandpass(1),theta_bandpass(2),8,0);
[~,mIdx]=max(wave);%get index max power for each timepiont
pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
power = rms(abs(wave))';

if use_ratioThetaDelta
    [wave_delta,~,~,~,~,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,delta_bandpass(1),delta_bandpass(2),8,0);
    power_delta = rms(abs(wave_delta))';

    power = power./power_delta;
end

if uselog10Power
    power = log10(power);
end

% find high noise periods
if ~isempty(threshold_noise) && threshold_noise
    M = movstd(double(lfpT.data),1 * lfpT.samplingRate);
    sw_std.t = downsample(lfpT.timestamps,1250); 
    sw_std.data = zscore(downsample(M,1250)); 
    sw_std.ints = [sw_std.t sw_std.t+1];
            
    intervals_below_threshold = sw_std.ints(find(sw_std.data<threshold_noise),:);
    clean_intervals = ConsolidateIntervals(intervals_below_threshold);
    clean_samples = InIntervals(lfpT.timestamps, clean_intervals);
    intervals = clean_intervals;
else
    clean_samples = ones(size(lfpT.timestamps));
    intervals = [lfpT.timestamps(1) lfpT.timestamps(end)];
end

disp('finding intervals below power threshold...')
thresh = mean(power(clean_samples)) + std(power(clean_samples))*powerThresh;
minWidth = (samplingRate./theta_bandpass(2)) * 3; % set the minimum width to four cycles

below=find(power<thresh);
below_thresh = [];
if max(diff(diff(below))) == 0
    below_thresh = [below(1) below(end)];
elseif length(below)>0;
    ends=find(diff(below)~=1);
    ends(end+1)=length(below);
    ends=sort(ends);
    lengths=diff(ends);
    stops=below(ends)./samplingRate;
    starts=lengths./samplingRate;
    starts = [1; starts];
    below_thresh(:,2)=stops;
    below_thresh(:,1)=stops-starts;
else
    below_thresh=[];
end
% now merge interval sets from input and power threshold
intervals = SubtractIntervals(intervals,below_thresh);  % subtract out low power intervals

if ~isempty(powerThreshold_nonTheta)
    disp('finding intervals below power threshold...')
    non_theta_power = power;
    non_theta_power = non_theta_power - mean(non_theta_power);
    non_theta_power = non_theta_power * -1;

    thresh = mean(non_theta_power(clean_samples)) + std(non_theta_power(clean_samples))*powerThreshold_nonTheta;
    minWidth = (samplingRate./theta_bandpass(2)) * 3; % set the minimum width to four cycles
    
    below=find(non_theta_power<thresh);
    below_thresh = [];
    if max(diff(diff(below))) == 0
        below_thresh = [below(1) below(end)];
    elseif length(below)>0;
        ends=find(diff(below)~=1);
        ends(end+1)=length(below);
        ends=sort(ends);
        lengths=diff(ends);
        stops=below(ends)./samplingRate;
        starts=lengths./samplingRate;
        starts = [1; starts];
        below_thresh(:,2)=stops;
        below_thresh(:,1)=stops-starts;
    else
        below_thresh=[];
    end
    % now merge interval sets from input and power threshold
    intervals_nonTheta = SubtractIntervals(clean_intervals,below_thresh);  % subtract out low power intervals 
    
    try
        disp('Discarting ripples intervals...');
        ripples = rippleMasterDetector;
        intervals_nonTheta_nonRipples = SubtractIntervals(intervals_nonTheta,ripples.timestamps);  % subtract out low power intervals 
        intervals_nonTheta_nonRipples(minWidth < diff(intervals_nonTheta_nonRipples')) = [];

    catch
        disp('Discarting ripples was not possible');
    end

end

thetaEpochs.lfpphase = lfpphase;
thetaEpochs.samplingRate = samplingRate;
thetaEpochs.power = power;
thetaEpochs.timestamps = t;
thetaEpochs.intervals = intervals;
thetaEpochs.channel = channel;
thetaEpochs.intervals_nonTheta = intervals_nonTheta;
thetaEpochs.intervals_nonTheta_nonRipples = intervals_nonTheta_nonRipples;
thetaEpochs.powerThreshold_nonTheta = powerThreshold_nonTheta;
thetaEpochs.params.powerThresh = powerThresh;
thetaEpochs.params.theta_bandpass = theta_bandpass;
thetaEpochs.params.powerThreshold_nonTheta = powerThreshold_nonTheta;
thetaEpochs.params.delta_bandpass = delta_bandpass;
thetaEpochs.params.use_ratioThetaDelta = use_ratioThetaDelta;
thetaEpochs.params.threshold_noise = threshold_noise;
thetaEpochs.params.uselog10Power = uselog10Power;
[thetaEpochs.idx.idx,thetaEpochs.idx.timestamps] = bz_INTtoIDX({thetaEpochs.intervals},'sf',1);
[thetaEpochs.idx_nonTheta.idx,thetaEpochs.idx_nonTheta.timestamps] = ...
    bz_INTtoIDX({thetaEpochs.intervals_nonTheta},'sf',1);
[thetaEpochs.idx_nonTheta_nonRipples.idx,thetaEpochs.idx_nonTheta_nonRipples.timestamps] = ...
    bz_INTtoIDX({thetaEpochs.intervals_nonTheta_nonRipples},'sf',1);

% try separating RUN and REM, and Quiet from Sleep in nonTheta
try SleepState = SleepScoreMaster(pwd,'noPrompts',true);
    
    % separating RUN and REM
    thetaRun_times = intersect(SleepState.idx.timestamps(SleepState.idx.states == 1),...
        thetaEpochs.idx.timestamps(thetaEpochs.idx.idx)); % 1 is WAKE
    thetaEpochs.thetaRun.idx = zeros(size(1:length(SleepState.idx.timestamps)));
    thetaEpochs.thetaRun.timestamps = SleepState.idx.timestamps;
    thetaEpochs.thetaRun.idx(thetaRun_times) = 1;
    thetaEpochs.thetaRun.states = thetaEpochs.thetaRun.idx;
    thetaEpochs.thetaRun.statenames = {'theta'};
    temp = bz_IDXtoINT(thetaEpochs.thetaRun);
    thetaEpochs.thetaRun.ints = temp.thetastate;
    
    thetaEpochs.thetaREM.timestamps = SleepState.idx.timestamps;
    thetaEpochs.thetaREM.idx = double(SleepState.idx.states==5);
    thetaEpochs.thetaREM.ints = SleepState.ints.REMstate;
    
    % isolating QWake
    QWake_times = intersect(SleepState.idx.timestamps(SleepState.idx.states == 1),...
        thetaEpochs.idx_nonTheta.timestamps(thetaEpochs.idx_nonTheta.idx)); % 1 is WAKE
    thetaEpochs.QWake.idx = zeros(size(1:length(SleepState.idx.timestamps)));
    thetaEpochs.QWake.timestamps = SleepState.idx.timestamps;
    thetaEpochs.QWake.idx(QWake_times) = 1;
    thetaEpochs.QWake.states = thetaEpochs.QWake.idx;
    thetaEpochs.QWake.statenames = {'QWake'};
    temp = bz_IDXtoINT(thetaEpochs.QWake);
    thetaEpochs.QWake.ints = temp.QWakestate;

    QWake_nonRipples_times = intersect(SleepState.idx.timestamps(SleepState.idx.states == 1),...
        thetaEpochs.idx_nonTheta_nonRipples.timestamps(thetaEpochs.idx_nonTheta_nonRipples.idx)); % 1 is WAKE
    thetaEpochs.QWake_nonRipples.idx = zeros(size(1:length(SleepState.idx.timestamps)));
    thetaEpochs.QWake_nonRipples.timestamps = SleepState.idx.timestamps;
    thetaEpochs.QWake_nonRipples.idx(QWake_nonRipples_times) = 1;
    thetaEpochs.QWake_nonRipples.states = thetaEpochs.QWake_nonRipples.idx;
    thetaEpochs.QWake_nonRipples.statenames = {'QWake_nonRipples'};
    temp = bz_IDXtoINT(thetaEpochs.QWake_nonRipples);
    thetaEpochs.QWake_nonRipples.ints = temp.QWake_nonRipplesstate;

    % isolating NREM
    NREM_times = intersect(SleepState.idx.timestamps(SleepState.idx.states == 3),...
        thetaEpochs.idx_nonTheta.timestamps(thetaEpochs.idx_nonTheta.idx)); % 3 is NREM
    thetaEpochs.NREM.idx = zeros(size(1:length(SleepState.idx.timestamps)));
    thetaEpochs.NREM.timestamps = SleepState.idx.timestamps;
    thetaEpochs.NREM.idx(NREM_times) = 1;
    thetaEpochs.NREM.states = thetaEpochs.NREM.idx;
    thetaEpochs.NREM.statenames = {'NREM'};
    temp = bz_IDXtoINT(thetaEpochs.NREM);
    thetaEpochs.NREM.ints = temp.NREMstate;

    NREM_nonRipples_times = intersect(SleepState.idx.timestamps(SleepState.idx.states == 3),...
        thetaEpochs.idx_nonTheta_nonRipples.timestamps(thetaEpochs.idx_nonTheta_nonRipples.idx)); % 3 is NREM
    thetaEpochs.NREM_nonRipples.idx = zeros(size(1:length(SleepState.idx.timestamps)));
    thetaEpochs.NREM_nonRipples.timestamps = SleepState.idx.timestamps;
    thetaEpochs.NREM_nonRipples.idx(NREM_nonRipples_times) = 1;
    thetaEpochs.NREM_nonRipples.states = thetaEpochs.NREM_nonRipples.idx;
    thetaEpochs.NREM_nonRipples.statenames = {'NREM_nonRipples'};
    temp = bz_IDXtoINT(thetaEpochs.NREM_nonRipples);
    thetaEpochs.NREM_nonRipples.ints = temp.NREM_nonRipplesstate;
catch 
    warning('Separating Run and REM was not possible!');
end

if updateSleepStates
    load([basenameFromBasepath(pwd) '.SleepState.states.mat'])
<<<<<<< HEAD
%     keyboard;
%     SleepState.ints.
%     thetaEpochs.thetaRun.idx
%     SleepState.detectorinfo
%     SleepState.ints.WAKEtheta2 = thetaEpochs.thetaRun
=======
        
    SleepState.ints.WAKEtheta_ThDt = thetaEpochs.thetaRun.ints;
    SleepState.ints.REMtheta_ThDt = thetaEpochs.thetaREM.ints;
    SleepState.ints.QWake_ThDt = thetaEpochs.QWake.ints;
    SleepState.ints.QWake_noRipples_ThDt = thetaEpochs.QWake_nonRipples.ints;
    SleepState.ints.NREM_ThDt = thetaEpochs.NREM.ints;
    SleepState.ints.NREM_noRipples_ThDt = thetaEpochs.NREM_nonRipples.ints;

    save([basenameFromBasepath(pwd) '.SleepState.states.mat'],'SleepState');
>>>>>>> 2dbe996b29da75ef9d12a43067f79fe96ee6bf21
end

if saveMat
    disp('Saving...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.thetaEpochs.states.mat'],'thetaEpochs');
end

if plotting
    params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
    [S,t,f] = mtspecgramc_fast(single(lfpT.data),[2 1],params); S(S==0) = NaN;
    S = log10(S); % in Db
    %S_det= bsxfun(@minus,S,polyval(polyfit(f,nanmean(S,1),2),f)); % detrending
    S_det= detrend(S',2)';

    figure;
    subplot(3,3,[1 2])
    imagesc(t,f,S_det',[-1.5 1.5]);
    ylim([1.5 50]);
    set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');

    subplot(3,3,[4 5])
    t_theta = sum(diff(thetaEpochs.intervals')) * diff(t);
    imagesc([0 t_theta],f,S_det(InIntervals(t,thetaEpochs.intervals),:)',[-1.5 1.5]);
    ylim([1.5 50]);
    set(gca,'TickDir','out'); ylabel('Theta epochs [Freq, Hz]');
    colormap jet

    subplot(3,3,[7 8])
    t_non_theta = sum(diff(thetaEpochs.intervals_nonTheta')) * diff(t);
    imagesc([0 t_non_theta],f,S_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:)',[-1.5 1.5]);
    ylim([1.5 50]);
    set(gca,'TickDir','out'); ylabel('Non-Theta epochs [Freq, Hz]');
    colormap jet
    
    subplot(2,3,[3 6])
    plotFill(f,S_det,'color', [.8 .8 .8],'lineStyle', '--'); xlim([1 30]);
    plotFill(f,S_det(InIntervals(t,thetaEpochs.intervals),:),'color', [.8 .2 .2],'lineStyle', '-'); xlim([1 30]);
    plotFill(f,S_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:),'color', [.2 .2 .8],'lineStyle', '-'); xlim([1 30]);
    ax = axis;
    fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
    ylabel('Full recording [Freq, Hz]'); xlabel('Freq [Hz]');   
    
    if saveSummary
        mkdir('SummaryFigures'); % create folder
        saveas(gcf,'SummaryFigures\thetaEpochs.png');
    end
end


cd(prevBasepath);
end