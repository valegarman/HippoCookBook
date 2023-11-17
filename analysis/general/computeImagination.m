function [semanticWords] = computeImagination(varargin)
% [semanticWords] = computeSemanticWords(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% for imagination words
%
% <OPTIONALS>
% spikes        buzcode spikes structure, if not provided tries loadSpikes.
% basepath      By default pwd.
% binSize       In seconds, default, 0.005.
% winSize       In seconds, default, 0.6.
% winSizePlot   Interval in seconds, default [-.3 .3] 
% plotOpt       Default true.
% saveFig       Default true
% saveMat       Default true.
% force         Overwrite previos analysis, default false.
% excludeStimulationPeriods
%               Default, true
% excludeIntervals 
%               2xN intervals to exlude
% winIndex      Default, [-.01 .01];
% interp0       Default, true.
% useBrainRegions
%               Default, true.
% includeIntervals
%               intervals to include for analysis
% OUTPUTS
% averageCCG
%
% Pablo Abad 2023
% 1: Objects
% 2: Spatial
% 3: Abstract
% 4: Actions

% Parse options
p = inputParser;
addParameter(p,'spikes',[]);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'binSize',0.005,@isnumeric);
addParameter(p,'winSize',0.6,@isnumeric);
addParameter(p,'plotOpt',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'winSizePlot',[-.3 .3],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'includeIntervals',[0 Inf],@isnumeric);
addParameter(p,'winIndex',[-.01 .01],@isnumeric);
addParameter(p,'interp0',[-.01 .01],@isnumeric);
addParameter(p,'useBrainRegions',true,@islogical);
addParameter(p,'useDistinctShanks',true,@islogical);
addParameter(p,'useCellType',true,@islogical);
addParameter(p,'channel',[],@isnumeric);

parse(p, varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
winSizePlot = p.Results.winSizePlot;
plotOpt = p.Results.plotOpt;
saveFig = p.Results.saveFig;
saveMat = p.Results.saveMat;
force = p.Results.force;
excludeIntervals = p.Results.excludeIntervals;
includeIntervals = p.Results.includeIntervals;
winIndex = p.Results.winIndex;
interp0 = p.Results.interp0;
useBrainRegions = p.Results.useBrainRegions;
useDistinctShanks = p.Results.useDistinctShanks;
useCellType = p.Results.useCellType;
channel = p.Results.channel;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.semanticWords.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('semanticWords already computed! Loading file...');
    load(targetFile.name);
    return
end

pulses = getAnalogPulses;

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if ~isempty(dir('semantic_words_6.mat'))
    targetFile = dir('semantic_words_6.mat');
    load(targetFile.name);
end

indexes = cell2mat(words(:,2));
indexes_2 = cell2mat(words(:,3));

count_even = 1;
count_odd = 1;
for ii = 1:size(indexes,1)
    if rem(ii,2) == 0
        even_indexes(count_even) = ii;
        count_even = count_even + 1;
    elseif rem(ii,2) == 1
        odd_indexes(count_odd) = ii;
        count_odd = count_odd + 1;
    end
end

abstract_indexes = find(indexes==3)';
action_indexes = find(indexes==4)';



%% PRIMER

psth_primer = spikesPsth([pulses.timestampsOn{channel}(even_indexes)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_primer.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_primer.responsecurveSmooth,1)
    psth_primer.responseZ(ii,:) = psth_primer.responsecurveZSmooth(ii,:) - ...
        mean(psth_primer.responsecurveSmooth(ii,win_Z))./std(psth_primer.responsecurveSmooth(ii,win_Z));
end
psth_primer.peakResponse = nanmean(psth_primer.responsecurve(:,win),2);
psth_primer.peakResponseZ = nanmean(psth_primer.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_primer.responsecurveZSmooth(find(~isnan(psth_primer.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primer.peakResponseZ(find(~isnan(psth_primer.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primer.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


%% SOUND

psth_sound = spikesPsth([pulses.timestampsOn{channel}(odd_indexes)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_sound.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_sound.responsecurveSmooth,1)
    psth_sound.responseZ(ii,:) = psth_sound.responsecurveZSmooth(ii,:) - ...
        mean(psth_sound.responsecurveSmooth(ii,win_Z))./std(psth_sound.responsecurveSmooth(ii,win_Z));
end
psth_sound.peakResponse = nanmean(psth_sound.responsecurve(:,win),2);
psth_sound.peakResponseZ = nanmean(psth_sound.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_sound.responsecurveZSmooth(find(~isnan(psth_sound.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_sound.peakResponseZ(find(~isnan(psth_sound.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_sound.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar

%% PRIMER ABSTRACT
a = ismember(even_indexes,abstract_indexes);
psth_primerabstract = spikesPsth([pulses.timestampsOn{channel}(even_indexes(a))'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_primerabstract.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_primerabstract.responsecurveSmooth,1)
    psth_primerabstract.responseZ(ii,:) = psth_primerabstract.responsecurveZSmooth(ii,:) - ...
        mean(psth_primerabstract.responsecurveSmooth(ii,win_Z))./std(psth_primerabstract.responsecurveSmooth(ii,win_Z));
end
psth_primerabstract.peakResponse = nanmean(psth_primerabstract.responsecurve(:,win),2);
psth_primerabstract.peakResponseZ = nanmean(psth_primerabstract.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_primerabstract.responsecurveZSmooth(find(~isnan(psth_primerabstract.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primerabstract.peakResponseZ(find(~isnan(psth_primerabstract.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primerabstract.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


%% PRIMER ACTION
a = ismember(even_indexes,action_indexes);
psth_primeraction = spikesPsth([pulses.timestampsOn{channel}(action_indexes(a))'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_primeraction.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_primeraction.responsecurveSmooth,1)
    psth_primeraction.responseZ(ii,:) = psth_primeraction.responsecurveZSmooth(ii,:) - ...
        mean(psth_primeraction.responsecurveSmooth(ii,win_Z))./std(psth_primeraction.responsecurveSmooth(ii,win_Z));
end
psth_primeraction.peakResponse = nanmean(psth_primeraction.responsecurve(:,win),2);
psth_primeraction.peakResponseZ = nanmean(psth_primeraction.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_primeraction.responsecurveZSmooth(find(~isnan(psth_primeraction.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primeraction.peakResponseZ(find(~isnan(psth_primeraction.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primeraction.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


%% SOUND ABSTRACT
a = ismember(odd_indexes,abstract_indexes);
psth_soundabstract = spikesPsth([pulses.timestampsOn{channel}(odd_indexes(a))'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_soundabstract.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_soundabstract.responsecurveSmooth,1)
    psth_soundabstract.responseZ(ii,:) = psth_soundabstract.responsecurveZSmooth(ii,:) - ...
        mean(psth_soundabstract.responsecurveSmooth(ii,win_Z))./std(psth_soundabstract.responsecurveSmooth(ii,win_Z));
end
psth_soundabstract.peakResponse = nanmean(psth_soundabstract.responsecurve(:,win),2);
psth_soundabstract.peakResponseZ = nanmean(psth_soundabstract.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_soundabstract.responsecurveZSmooth(find(~isnan(psth_soundabstract.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_soundabstract.peakResponseZ(find(~isnan(psth_soundabstract.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_soundabstract.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar

%% SOUND ACTION
a = ismember(odd_indexes,action_indexes);
psth_soundaction = spikesPsth([pulses.timestampsOn{channel}(odd_indexes(a))'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_soundaction.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_soundaction.responsecurveSmooth,1)
    psth_soundaction.responseZ(ii,:) = psth_soundaction.responsecurveZSmooth(ii,:) - ...
        mean(psth_soundaction.responsecurveSmooth(ii,win_Z))./std(psth_soundaction.responsecurveSmooth(ii,win_Z));
end
psth_soundaction.peakResponse = nanmean(psth_soundaction.responsecurve(:,win),2);
psth_soundaction.peakResponseZ = nanmean(psth_soundaction.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_soundaction.responsecurveZSmooth(find(~isnan(psth_soundaction.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_soundaction.peakResponseZ(find(~isnan(psth_soundaction.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_soundaction.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar

%% COHERENT

psth_coherent = spikesPsth([pulses.timestampsOn{channel}(indexes_2 == 1)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_primer.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_primer.responsecurveSmooth,1)
    psth_primer.responseZ(ii,:) = psth_primer.responsecurveZSmooth(ii,:) - ...
        mean(psth_primer.responsecurveSmooth(ii,win_Z))./std(psth_primer.responsecurveSmooth(ii,win_Z));
end
psth_primer.peakResponse = nanmean(psth_primer.responsecurve(:,win),2);
psth_primer.peakResponseZ = nanmean(psth_primer.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_primer.responsecurveZSmooth(find(~isnan(psth_primer.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primer.peakResponseZ(find(~isnan(psth_primer.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primer.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


%% NONCOHERENT

psth_primer = spikesPsth([pulses.timestampsOn{channel}(even_indexes)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);
win_resp = [0.8 1.1];
ts_psth = psth_primer.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_primer.responsecurveSmooth,1)
    psth_primer.responseZ(ii,:) = psth_primer.responsecurveZSmooth(ii,:) - ...
        mean(psth_primer.responsecurveSmooth(ii,win_Z))./std(psth_primer.responsecurveSmooth(ii,win_Z));
end
psth_primer.peakResponse = nanmean(psth_primer.responsecurve(:,win),2);
psth_primer.peakResponseZ = nanmean(psth_primer.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_primer.responsecurveZSmooth(find(~isnan(psth_primer.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primer.peakResponseZ(find(~isnan(psth_primer.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primer.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar
end