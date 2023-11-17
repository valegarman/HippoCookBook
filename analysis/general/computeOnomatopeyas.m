function [onomatopeyas] = computeOnomatopeyas(varargin)
% [semanticWords] = computeSemanticWords(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% for semantic words
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

targetFile = dir('*.onomatopeyas.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('semanticWords already computed! Loading file...');
    load(targetFile.name);
    return
end

pulses = getAnalogPulses;

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if ~isempty(dir('onomatopeyas.mat'))
    targetFile = dir('onomatopeyas.mat');
    load(targetFile.name);
end

indexes = cell2mat(onomatopeyas(:,2:3));

%% ============== PRIMER VS STIMULI ===========

% Primer
psth_primer = spikesPsth([pulses.timestampsOn{channel}(indexes(:,1) == 1)'],'numRep',100,'saveMat',false,...
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

% figure;
% imagesc_ranked(ts_psth,[],psth_primer.responsecurveZSmooth(find(~isnan(psth_primer.responsecurveZSmooth(:,1))),:),[-3 3],...
%     psth_primer.peakResponseZ(find(~isnan(psth_primer.responsecurveZSmooth(:,1)))));
% colormap jet;
% hold on;
% zmean = nanmean(zscore(psth_primer.responsecurveZSmooth,[],2));
% zmean = zmean - min(zmean); 
% zmean = zmean/max(zmean) * (111-1) * std(zmean);
% plot(ts_psth,zmean+10,'k','LineWidth',1);
% colorbar


figure;
imagesc_ranked(ts_psth,[],psth_primer.responsecurveZSmooth(find(~isnan(psth_primer.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_primer.peakResponseZ(find(~isnan(psth_primer.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_primer.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_primer = (length(find(psth_primer.bootsTrapTest == 1)) / length(psth_primer.bootsTrapTest))*100;
negative_responsive_cells_primer = (length(find(psth_primer.bootsTrapTest == -1)) / length(psth_primer.bootsTrapTest))*100;

% Stimuli
psth_stim = spikesPsth(pulses.timestampsOn{channel}(indexes(:,1) == 2)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_stim.responsecurveSmooth,1)
    psth_stim.responseZ(ii,:) = psth_stim.responsecurveZSmooth(ii,:) - ...
        mean(psth_stim.responsecurveSmooth(ii,win_Z))./std(psth_stim.responsecurveSmooth(ii,win_Z));
end
psth_stim.peakResponse = nanmean(psth_stim.responsecurve(:,win),2);
psth_stim.peakResponseZ = nanmean(psth_stim.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_stim.responsecurveZSmooth(find(~isnan(psth_stim.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_stim.peakResponseZ(find(~isnan(psth_stim.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_stim.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_stim = (length(find(psth_stim.bootsTrapTest == 1)) / length(psth_stim.bootsTrapTest))*100;
negative_responsive_cells_stim = (length(find(psth_stim.bootsTrapTest == -1)) / length(psth_stim.bootsTrapTest))*100;



%% ======= COHERENTE VS NON-COHERENT ============

% Coherent
psth_coherent = spikesPsth([pulses.timestampsOn{channel}(indexes(:,1) == 2 & indexes(:,2) == 1)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

win_resp = [0.8 1.1];
ts_psth = psth_coherent.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_coherent.responsecurveSmooth,1)
    psth_coherent.responseZ(ii,:) = psth_coherent.responsecurveZSmooth(ii,:) - ...
        mean(psth_coherent.responsecurveSmooth(ii,win_Z))./std(psth_coherent.responsecurveSmooth(ii,win_Z));
end
psth_coherent.peakResponse = nanmean(psth_coherent.responsecurve(:,win),2);
psth_coherent.peakResponseZ = nanmean(psth_coherent.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_coherent.responsecurveZSmooth(find(~isnan(psth_coherent.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_coherent.peakResponseZ(find(~isnan(psth_coherent.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_coherent.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


figure;
imagesc_ranked(ts_psth,[],psth_coherent.responsecurveZSmooth(find(~isnan(psth_coherent.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_coherent.peakResponseZ(find(~isnan(psth_coherent.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_coherent.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_coherent = (length(find(psth_coherent.bootsTrapTest == 1)) / length(psth_coherent.bootsTrapTest))*100;
negative_responsive_cells_coherent = (length(find(psth_coherent.bootsTrapTest == -1)) / length(psth_coherent.bootsTrapTest))*100;

% NonCoherent
psth_noncoherent = spikesPsth(pulses.timestampsOn{channel}(indexes(:,1) == 2 & indexes(:,2) == 0)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_noncoherent.responsecurveSmooth,1)
    psth_noncoherent.responseZ(ii,:) = psth_noncoherent.responsecurveZSmooth(ii,:) - ...
        mean(psth_noncoherent.responsecurveSmooth(ii,win_Z))./std(psth_noncoherent.responsecurveSmooth(ii,win_Z));
end
psth_noncoherent.peakResponse = nanmean(psth_noncoherent.responsecurve(:,win),2);
psth_noncoherent.peakResponseZ = nanmean(psth_noncoherent.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_noncoherent.responsecurveZSmooth(find(~isnan(psth_noncoherent.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_noncoherent.peakResponseZ(find(~isnan(psth_noncoherent.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_noncoherent.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_noncoherent = (length(find(psth_noncoherent.bootsTrapTest == 1)) / length(psth_noncoherent.bootsTrapTest))*100;
negative_responsive_cells_noncoherent = (length(find(psth_noncoherent.bootsTrapTest == -1)) / length(psth_noncoherent.bootsTrapTest))*100;




onomatopeyas = [];





end