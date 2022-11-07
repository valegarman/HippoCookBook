function [cohgram] = computeCohgram(varargin)
% Computes coherogram between two channels.
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
% 'bandpass'            Default []
% 'channel'             Numeric [ex, 5]; by default calls
%                           brain Regions.
% 'updateSleepStates'   Default true
% 'useCSD'              Default, true.
% 'discardRipples'      Discard ripples from nonTheta, default true.
% 'useSubfolders'       Computes PowerSpectrum for each subfolder
% 'useThetaEpochs'      Compute PowerSpectrum only for detected Theta
%                       Epochs
% 
% OUTPUT
% 'cohgram'               structure with power spectrum and coherence between
%                       channels
%
% Pablo Abad 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO
% 1. Include options for more channels
% ...
%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp1',[],@isstruct);
addParameter(p,'lfp2',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'passband',[0 200], @isnumeric);
addParameter(p,'powerThreshold',[.8], @isnumeric);
addParameter(p,'channel1',[],@isnumeric);
addParameter(p,'channel2',[],@isnumeric);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'updateSleepStates',true,@islogical);
addParameter(p,'useCSD',false,@islogical);
addParameter(p,'delta_bandpass',[1 3],@isnumeric);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'lgamma_bandpass',[20 60],@isnumeric);
addParameter(p,'hgamma_bandpass',[60 100],@isnumeric);
addParameter(p,'use_ratioThetaDelta',true,@islogical);
addParameter(p,'threshold_noise',2.5,@isnumeric);
addParameter(p,'uselog10Power',true,@islogical);
addParameter(p,'powerThreshold_nonTheta',.5,@isnumeric);
addParameter(p,'discardRipples',true,@islogical);
addParameter(p,'useSubFolders',true,@islogical);
addParameter(p,'useThetaEpochs',true,@islogical);
addParameter(p,'useRippleChannel',true,@islogical);
addParameter(p,'useThetaChannel',false,@islogical);
addParameter(p,'notchFilter',false,@islogical);
addParameter(p,'plt',true,@islogical);
addParameter(p,'saveFig',true,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp1 = p.Results.lfp1;
lfp2 = p.Results.lfp2;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
passband = p.Results.passband;
powerThresh = p.Results.powerThreshold;
channel1 = p.Results.channel1;
channel2 = p.Results.channel2;
plotting = p.Results.plotting;
updateSleepStates = p.Results.updateSleepStates;
useCSD = p.Results.useCSD;
delta_bandpass = p.Results.delta_bandpass;
theta_bandpass = p.Results.theta_bandpass;
lgamma_bandpass = p.Results.lgamma_bandpass;
hgamma_bandpass = p.Results.hgamma_bandpass;
use_ratioThetaDelta = p.Results.use_ratioThetaDelta;
threshold_noise =p.Results.threshold_noise;
uselog10Power = p.Results.uselog10Power;
powerThreshold_nonTheta = p.Results.powerThreshold_nonTheta;
discardRipples = p.Results.discardRipples;
useSubFolders = p.Results.useSubFolders;
useThetaEpochs = p.Results.useThetaEpochs;
useRippleChannel = p.Results.useRippleChannel;
useThetaChannel = p.Results.useThetaChannel;
notchFilter = p.Results.notchFilter;
plt = p.Results.plt;
saveFig = p.Results.saveFig;

% Deal with inputs
prevBasepath = pwd;
cd(basepath);

session = loadSession(basepath);

targetFile = dir('*.cohgram.lfp.mat');
if ~isempty(targetFile) && ~force
    disp('Power Spectrum already detected! Loading file.');
    load(targetFile.name);
    return
end

targetFile = dir('*.ripples.events.mat');
if ~isempty(targetFile)
    disp('Ripples already detected. loading file...');
    load(targetFile.name);
end

% Load brainRegions to load channels
targetFile = dir('*.brainRegions.channelinfo.mat');
if ~isempty(targetFile)
    disp('Brain Regions already detected ! Loading file...');
    load(targetFile.name);
    regions = fields(brainRegions);
end

% Load MergePoints
if ~isempty(dir([session.general.name,'.MergePoints.events.mat']))
    disp('MergePoints detected. Loading file...');
    file = dir([session.general.name,'.MergePoints.events.mat']);
    load(file.name);
end

% Load tracking
if ~isempty(dir([session.general.name,'.Tracking.Behavior.mat']))
    disp('Tracking detected. Loading file');
    file = dir([session.general.name,'.Tracking.Behavior.mat']);
    load(file.name);
end

% Load theta Epochs
if ~isempty(dir([session.general.name,'.thetaEpochs.states.mat']))
    disp('Theta Epochs detected. Loading file');
    file = dir([session.general.name,'.thetaEpochs.states.mat']);
    load(file.name);
end
    
% Remove from brain regions channels that are defined as bad

for ii = 1:length(regions)
    brainRegions.(regions{ii}).channels(ismember(brainRegions.(regions{ii}).channels,session.channelTags.Bad.channels)) = [];
end

% Get as many channels as brain regions are defined
if isempty(channel1) && isempty(channel2)
    for ii = 1:length(regions)
        channel(ii) = brainRegions.(regions{ii}).channels(randperm(length(brainRegions.(regions{ii}).channels),1));
        lfp(ii,:) = getLFP(channel(ii),'noPrompts',true);     
    end
end

params.Fs = session.extracellular.srLfp; params.fpass = [1 200]; params.tapers = [3 5]; params.pad = 1;
channels_pairs = nchoosek(1:length(lfp), 2);

% Computes cohgram for whole recording
for ii = 1:size(channels_pairs,1)
    if notchFilter
        lfp1 = Notched(single(lfp(channels_pairs(ii,1)).data),params.Fs,50,100,150,200);
        lfp2 = Notched(single(lfp(channels_pairs(ii,2)).data),params.Fs,50,100,150,200);
    else
        lfp1 = single(lfp(channels_pairs(ii,1)).data);
        lfp2 = single(lfp(channels_pairs(ii,2)).data);
    end
    [coherogram,phase,S12,S1,S2,t,f] = cohgramc(lfp1,lfp2,[2 1], params);
    
    coherogram(coherogram==0) = NaN;
    phase(phase==0) = NaN;
    S12(S12==0) = NaN;
    S1(S1==0) = NaN;
    S2(S2==0) = NaN;
    
    S12 = log10(S12); % in Db
    S1 = log10(S1); % in Db
    S2 = log10(S2);
    
    S12_det = detrend(S12',2)';
    S1_det = detrend(S1',2)';
    S2_det = detrend(S2',2)';
    coherogram = detrend(coherogram',2)';
    phase = detrend(phase',2)';
    
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).coherogram = coherogram;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).phase = phase;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).t = t;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).f = f;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).S1 = S1_det;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).S2 = S2_det;
    
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).ch1 = lfp(channels_pairs(ii,1)).channels;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).ch2 = lfp(channels_pairs(ii,2)).channels;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).region1 = lfp(channels_pairs(ii,1)).region;
    cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).region2 = lfp(channels_pairs(ii,2)).region;
    
    if useThetaEpochs
        % Theta Epochs
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.coherogram = coherogram(InIntervals(t,thetaEpochs.intervals),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.phase = phase(InIntervals(t,thetaEpochs.intervals),:);
%         cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.t = t(InIntervals(t,thetaEpochs.intervals),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.S1 = S1_det(InIntervals(t,thetaEpochs.intervals),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.S2 = S2_det(InIntervals(t,thetaEpochs.intervals),:);
        
        % Non Theta Epochs
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).NonthetaEpochs.coherogram = coherogram(InIntervals(t,thetaEpochs.intervals_nonTheta),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).NonthetaEpochs.phase = phase(InIntervals(t,thetaEpochs.intervals_nonTheta),:);
%         cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).thetaEpochs.t = t(InIntervals(t,thetaEpochs.intervals),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).NonthetaEpochs.S1 = S1_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:);
        cohgram.([lfp(channels_pairs(ii,1)).region,'_',lfp(channels_pairs(ii,2)).region]).NonthetaEpochs.S2 = S2_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:);
    end
end

%% Plotting
if plt
    flds = fields(cohgram);
    for ii = 1:length(flds)
        t = cohgram.(flds{ii}).t;
        f = cohgram.(flds{ii}).f;
        % Figure 1 (Whole Recording)
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(4,6,[1 2])
        imagesc(t,f,cohgram.(flds{ii}).coherogram',[-1 1]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);  
        
        subplot(4,6,[7 8])
        imagesc(t,f,cohgram.(flds{ii}).phase',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);    
        
        subplot(4,6,[13 14])
        imagesc(t,f,cohgram.(flds{ii}).S1',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);  
        
        subplot(4,6,[19 20])
        imagesc(t,f,cohgram.(flds{ii}).S2',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        subplot(4,6,[3 9 15 21])
        plotFill(f,mean(cohgram.(flds{ii}).coherogram),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Full recording [r]'); xlabel('Freq [Hz]');  
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[4 10 16 22])
        plotFill(f,mean(cohgram.(flds{ii}).phase),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Full recording [Phase]'); xlabel('Freq [Hz]');  
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[5 11 17 23])
        plotFill(f,mean(cohgram.(flds{ii}).S1),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Full recording [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);
        
        subplot(4,6,[6 12 18 24])
        plotFill(f,mean(cohgram.(flds{ii}).S2),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Full recording [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        if saveFig
            saveas(gcf,['SummaryFigures\cohgram_',num2str(cohgram.(flds{ii}).region1),'_',num2str(cohgram.(flds{ii}).region2),'.png']);
        end
        
        
        % Figure 2 (Theta epochs)
        t_theta = sum(diff(thetaEpochs.intervals')) * diff(t);
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(4,6,[1 2])
        imagesc([0 t_theta],f,cohgram.(flds{ii}).thetaEpochs.coherogram',[-1 1]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[7 8])
        imagesc([0 t_theta],f,cohgram.(flds{ii}).thetaEpochs.phase',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);    
        
        subplot(4,6,[13 14])
        imagesc([0 t_theta],f,cohgram.(flds{ii}).thetaEpochs.S1',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);  
        
        subplot(4,6,[19 20])
        imagesc([0 t_theta],f,cohgram.(flds{ii}).thetaEpochs.S2',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        subplot(4,6,[3 9 15 21])
        plotFill(f,mean(cohgram.(flds{ii}).thetaEpochs.coherogram),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Theta Epochs [r]'); xlabel('Freq [Hz]');  
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[4 10 16 22])
        plotFill(f,mean(cohgram.(flds{ii}).thetaEpochs.phase),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Theta Epochs [Phase]'); xlabel('Freq [Hz]');  
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[5 11 17 23])
        plotFill(f,mean(cohgram.(flds{ii}).thetaEpochs.S1),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Theta Epochs [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);
        
        subplot(4,6,[6 12 18 24])
        plotFill(f,mean(cohgram.(flds{ii}).thetaEpochs.S2),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Theta Epochs [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        if saveFig
            saveas(gcf,['SummaryFigures\cohgram_thetaEpochs_',num2str(cohgram.(flds{ii}).region1),'_',num2str(cohgram.(flds{ii}).region2),'.png']);
        end
        
        % Figure 3 (Non-theta epochs)
        
        t_non_theta = sum(diff(thetaEpochs.intervals_nonTheta')) * diff(t);
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(4,6,[1 2])
        imagesc([0 t_non_theta],f,cohgram.(flds{ii}).NonthetaEpochs.coherogram',[-1 1]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[7 8])
        imagesc([0 t_non_theta],f,cohgram.(flds{ii}).NonthetaEpochs.phase',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);    
        
        subplot(4,6,[13 14])
        imagesc([0 t_non_theta],f,cohgram.(flds{ii}).NonthetaEpochs.S1',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);  
        
        subplot(4,6,[19 20])
        imagesc([0 t_non_theta],f,cohgram.(flds{ii}).NonthetaEpochs.S2',[-1.5 1.5]);
        colormap jet
        set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        subplot(4,6,[3 9 15 21])
        plotFill(f,mean(cohgram.(flds{ii}).NonthetaEpochs.coherogram),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Non Theta Epochs [r]'); xlabel('Freq [Hz]');  
        title(['Coherence (r) Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[4 10 16 22])
        plotFill(f,mean(cohgram.(flds{ii}).NonthetaEpochs.phase),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Non Theta Epochs [Phase]'); xlabel('Freq [Hz]');  
        title(['Phase Coherence Ch: ', num2str(cohgram.(flds{ii}).ch1) , ' Ch: ', num2str(cohgram.(flds{ii}).ch2)]);
        
        subplot(4,6,[5 11 17 23])
        plotFill(f,mean(cohgram.(flds{ii}).NonthetaEpochs.S1),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Non Theta Epochs [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch1), ' Region: ', cohgram.(flds{ii}).region1]);
        
        subplot(4,6,[6 12 18 24])
        plotFill(f,mean(cohgram.(flds{ii}).NonthetaEpochs.S2),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);
        ax = axis;
        fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
        fill([lgamma_bandpass flip(lgamma_bandpass)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
        fill([hgamma_bandpass flip(hgamma_bandpass)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
        ylabel('Non Theta Epochs [Freq, Hz]'); xlabel('Freq [Hz]');  
        title(['Ch: ', num2str(cohgram.(flds{ii}).ch2), ' Region: ', cohgram.(flds{ii}).region2]);
        
        if saveFig
            saveas(gcf,['SummaryFigures\cohgram_NonthetaEpochs_',num2str(cohgram.(flds{ii}).region1),'_',num2str(cohgram.(flds{ii}).region2),'.png']);
        end
        
    end
end

if saveMat
    try
        save([session.general.name,'.cohgram.lfp.mat','cohgram']);
    catch
        disp('Saving with -v7.3...');
        save([session.general.name,'.cohgram.lfp.mat','cohgram'],'-v7.3');
    end
end

cd(prevBasepath);

end

