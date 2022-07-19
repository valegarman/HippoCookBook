function [powerSpectrum] = computePowerSpectrum(varargin)
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
% 'bandpass'            Default []
% 'channel'             Numeric [ex, 5]; by default calls
%                           getHippocampalLayers and uses oriens.
% 'updateSleepStates'   Default true
% 'useCSD'              Default, true.
% 'discardRipples'      Discard ripples from nonTheta, default true.
% 'useSubfolders'       Computes PowerSpectrum for each subfolder
% 'useThetaEpochs'      Compute PowerSpectrum only for detected Theta
%                       Epochs
% 
% OUTPUT
% powerSpectrum        states structure with theta epochs intervals
%
% Pablo Abad 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'passband',[0 200], @isnumeric);
addParameter(p,'powerThreshold',[.8], @isnumeric);
addParameter(p,'channel',[],@isnumeric);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'updateSleepStates',true,@islogical);
addParameter(p,'useCSD',false,@islogical);
addParameter(p,'delta_bandpass',[1 3],@isnumeric);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'use_ratioThetaDelta',true,@islogical);
addParameter(p,'threshold_noise',2.5,@isnumeric);
addParameter(p,'uselog10Power',true,@islogical);
addParameter(p,'powerThreshold_nonTheta',.5,@isnumeric);
addParameter(p,'discardRipples',true,@islogical);
addParameter(p,'useSubFolders',true,@islogical);
addParameter(p,'useThetaEpochs',true,@islogical);
addParameter(p,'useRippleChannel',true,@islogical);
addParameter(p,'useThetaChannel',false,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
passband = p.Results.passband;
powerThresh = p.Results.powerThreshold;
channel = p.Results.channel;
plotting = p.Results.plotting;
updateSleepStates = p.Results.updateSleepStates;
useCSD = p.Results.useCSD;
delta_bandpass = p.Results.delta_bandpass;
theta_bandpass = p.Results.theta_bandpass;
use_ratioThetaDelta = p.Results.use_ratioThetaDelta;
threshold_noise =p.Results.threshold_noise;
uselog10Power = p.Results.uselog10Power;
powerThreshold_nonTheta = p.Results.powerThreshold_nonTheta;
discardRipples = p.Results.discardRipples;
useSubFolders = p.Results.useSubFolders;
useThetaEpochs = p.Results.useThetaEpochs;
useRippleChannel = p.Results.useRippleChannel;
useThetaChannel = p.Results.useThetaChannel;


% Deal with inputs
prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.powerSpectrum.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Theta epochs already detected! Loading file.');
    load(targetFile.name);
    return
end

targetFile = dir('*.ripples.events.mat');
if ~isempty(targetFile)
    disp('Ripples already detected. loading file...');
    load(targetFile.name);
end

session = loadSession(basepath);

% 
if isempty(channel)
    hippocampalLayers = getHippocampalLayers;
    if useRippleChannel
        % Get ripples
        channel = ripples.detectorinfo.detectionchannel;
    elseif useThetaChannel
        % first, try to load slm channel in best shank
        if ~isnan(hippocampalLayers.bestShankLayers.slm)
            channel = hippocampalLayers.bestShankLayers.slm;
        else
        % then, try to load it from any shank, and get the slm channel with the
        % highest theta power
        powerProfile_theta = powerSpectrumProfile([6 12]);
            for ii = 1:length(hippocampalLayers.layers)
                temp(ii) = hippocampalLayers.layers{ii}.slm;
            end
            if any(~isnan(temp))
                temp(isnan(temp)) = [];
                [~,idx]= max(powerProfile_theta.mean(ismember(powerProfile_theta.channels, temp)));
                channel = temp(idx);
            else
                % if not possible, just take channel with most theta power
                [~,idx]= max(powerProfile_theta.mean);
                channel = powerProfile_theta.channels(idx);
            end
        end
    end
end

if isempty(lfp) && ~useCSD
    lfp = getLFP(channel,'noPrompts',true);
elseif useCSD
    disp('Computing CSD...');
    lfp = computeCSD(lfp,'channels',channel);
else
    warning('CSD estimation not possible. Using LFP...');
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
    

samplingRate = lfp.samplingRate;
params.Fs = lfp.samplingRate; params.fpass = [1 200]; params.tapers = [3 5]; params.pad = 1;
[S,t,f] = mtspecgramc_fast(single(lfp.data), [2 1], params); S(S == 0) = NaN;
S = log10(S); % in Db
S_det = detrend(S',2)';


figure,
subplot(3,3,[1 2])
imagesc(t,f,S_det',[-1.5 1.5]);
ylim([1.5 200]);
set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
colormap jet

subplot(3,3,[4 5])
try
    t_theta = sum(diff(thetaEpochs.intervals')) * diff(t);
    imagesc([0 t_theta],f,S_det(InIntervals(t,thetaEpochs.intervals),:)',[-1.5 1.5]);
    ylim([1.5 200]);
    set(gca,'TickDir','out'); ylabel('Theta epochs [Freq, Hz]');
    colormap jet
end
    
subplot(3,3,[7 8])
try
    t_non_theta = sum(diff(thetaEpochs.intervals_nonTheta')) * diff(t);
    imagesc([0 t_non_theta],f,S_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:)',[-1.5 1.5]);
    ylim([1.5 200]);
    set(gca,'TickDir','out'); ylabel('Non-Theta epochs [Freq, Hz]');
    colormap jet
end

subplot(2,3,[3 6])
try
    plotFill(f,S_det,'color', [.8 .8 .8],'lineStyle', '--'); xlim([1 200]);
    plotFill(f,S_det(InIntervals(t,thetaEpochs.intervals),:),'color', [.8 .2 .2],'lineStyle', '-'); xlim([1 200]);
    plotFill(f,S_det(InIntervals(t,thetaEpochs.intervals_nonTheta),:),'color', [.2 .2 .8],'lineStyle', '-'); xlim([1 200]);
    ax = axis;
    fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
    ylabel('Full recording [Freq, Hz]'); xlabel('Freq [Hz]');   
end 
% if saveSummary
%     mkdir('SummaryFigures'); % create folder
%     saveas(gcf,'SummaryFigures\PowerSpectrum_WholeRecording.png');
% end 

%% Computing for Baseline vs Drug condition

% Open Field
openField_index = [];
tracking_index = find(ismember(MergePoints.foldernames,tracking.folders));
count = 1;
for ii = 1:length(tracking.folders)
    if strcmpi(tracking.apparatus{ii}.name,'Open Field')
        openField_index = [openField_index ii];
        t_openField{count} = MergePoints.timestamps(tracking_index(ii),:);
        count = count + 1;
    end
end

t_openField1 = sum(diff(t_openField{1}')) * diff(t);
t_openField2 = sum(diff(t_openField{2}')) * diff(t);


figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,3,[1 2])
imagesc([0 t_openField1],f,S_det(InIntervals(t,t_openField{1}),:)',[-1.5 1.5]);
ylim([1.5 200]);
set(gca,'TickDir','out'); ylabel('Open Field 1 [Freq, Hz]');
colormap jet, colorbar

subplot(3,3,[4 5])
imagesc([0 t_openField2],f,S_det(InIntervals(t,t_openField{2}),:)',[-1.5 1.5]);
ylim([1.5 200]);
set(gca,'TickDir','out'); ylabel('Open Field 2 [Freq, Hz]');
colormap jet, colorbar

subplot(3,3,7)
plot(tracking.position.x(InIntervals(tracking.timestamps,t_openField{1})),tracking.position.y(InIntervals(tracking.timestamps,t_openField{1})));
axis ij
title('Open Field 1');
xlabel('norm/cm'); ylabel('norm/cm');

subplot(3,3,8)
plot(tracking.position.x(InIntervals(tracking.timestamps,t_openField{2})),tracking.position.y(InIntervals(tracking.timestamps,t_openField{2})));
axis ij
title('Open Field 2');
xlabel('norm/cm'); ylabel('norm/cm');

subplot(3,3,[3 6 9])
plotFill(f,S_det(InIntervals(t,t_openField{1}),:),'color',[.8 .8 .8], 'lineStyle','-'); xlim([1 200]);
plotFill(f,S_det(InIntervals(t,t_openField{2}),:),'color',[1 0 0], 'lineStyle','-'); xlim([1 200]);
ax = axis;
fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
ylabel('OpenField1 - OpenField2 [Freq, Hz]'); xlabel('Freq [Hz]');

if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,'SummaryFigures\PowerSpectrum_OpenField.png');
end


% YMaze
YMaze_index = [];
count = 1;
for ii = 1:length(tracking.folders)
    if strcmpi(tracking.apparatus{ii}.name,'YMaze Apparatus')
        YMaze_index = [YMaze_index ii];
        t_YMaze{count} = MergePoints.timestamps(tracking_index(ii),:);
        count = count + 1;
    end
end

t_YMaze1 = sum(diff(t_YMaze{1}')) * diff(t);
t_YMaze2 = sum(diff(t_YMaze{2}')) * diff(t);


figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,3,[1 2])
imagesc([0 t_YMaze1],f,S_det(InIntervals(t,t_YMaze{1}),:)',[-1.5 1.5]);
ylim([1.5 200]);
set(gca,'TickDir','out'); ylabel('YMaze 1 [Freq, Hz]');
colormap jet, colorbar

subplot(3,3,[4 5])
imagesc([0 t_YMaze2],f,S_det(InIntervals(t,t_YMaze{2}),:)',[-1.5 1.5]);
ylim([1.5 200]);
set(gca,'TickDir','out'); ylabel('YMaze 2 [Freq, Hz]');
colormap jet, colorbar

subplot(3,3,7)
plot(tracking.position.x(InIntervals(tracking.timestamps,t_YMaze{1})),tracking.position.y(InIntervals(tracking.timestamps,t_YMaze{1})));
axis ij
title('YMaze 1');
xlabel('norm/cm'); ylabel('norm/cm');

subplot(3,3,8)
plot(tracking.position.x(InIntervals(tracking.timestamps,t_YMaze{2})),tracking.position.y(InIntervals(tracking.timestamps,t_YMaze{2})));
axis ij
title('YMaze 2');
xlabel('norm/cm'); ylabel('norm/cm');

subplot(3,3,[3 6 9])
plotFill(f,S_det(InIntervals(t,t_YMaze{1}),:),'color',[.8 .8 .8], 'lineStyle','-'); xlim([1 200]);
plotFill(f,S_det(InIntervals(t,t_YMaze{2}),:),'color',[1 0 0], 'lineStyle','-'); xlim([1 200]);
ax = axis;
fill([theta_bandpass flip(theta_bandpass)],[ax([3 3 4 4])],[.8 .5 .5],'EdgeColor','none','FaceAlpha',.1);
ylabel('YMaze1 - YMaze2 [Freq, Hz]'); xlabel('Freq [Hz]');

if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,'SummaryFigures\PowerSpectrum_YMaze.png');
end

close all;



%% OUTPUT
powerSpectrum = [];
% Subfolders analysis
if useSubFolders
    if ~isempty(dir([session.general.name,'.MergePoints.events.mat']))
        disp('MergePoints detected. Loading file...');
        file = dir([session.general.name,'.MergePoints.events.mat']);
        load(file.name);
    end
    
    for ii = 1:length(MergePoints.foldernames)
        powerSpectrum.(MergePoints.foldernames{ii}).S = S_det(InIntervals(t,MergePoints.timestamps(ii,:)));
        powerSpectrum.(MergePoints.foldernames{ii}).t = t(InIntervals(t,MergePoints.timestamps(ii,:)));
        powerSpectrum.(MergePoints.foldernames{ii}).f = f;
    end
end

powerSpectrum.t = t;
powerSpectrum.f = f;
powerSpectrum.S = S_det;
powerSpectrum.processinginfo.function = 'computePowerSpectrum';
powerSpectrum.processinginfo.params.frange = passband;


if saveMat
    disp('Saving...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.powerSpectrum.states.mat'],'powerSpectrum');
end

cd(prevBasepath);
end