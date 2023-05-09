function [semanticWords] = computeSemanticWords(varargin)
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

if ~isempty(dir('semantic_words.mat'))
    targetFile = dir('semantic_words.mat');
    load(targetFile.name);
end

indexes = cell2mat(words(:,2));

%% ============== GENERAL ALL WORDS ===========

psth_general = spikesPsth([pulses.timestampsOn{3}'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

win_resp = [0.8 1.1];
ts_psth = psth_general.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_general.responsecurveSmooth,1)
    psth_general.responseZ(ii,:) = psth_general.responsecurveZSmooth(ii,:) - ...
        mean(psth_general.responsecurveSmooth(ii,win_Z))./std(psth_general.responsecurveSmooth(ii,win_Z));
end
psth_general.peakResponse = nanmean(psth_general.responsecurve(:,win),2);
psth_general.peakResponseZ = nanmean(psth_general.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_general.responsecurveZSmooth(find(~isnan(psth_general.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_general.peakResponseZ(find(~isnan(psth_general.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_general.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);
colorbar


figure;
imagesc_ranked(ts_psth,[],psth_general.responsecurveZSmooth(find(~isnan(psth_general.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_general.peakResponseZ(find(~isnan(psth_general.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_general.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_general = (length(find(psth_general.bootsTrapTest == 1)) / length(psth_general.bootsTrapTest))*100;
negative_responsive_cells_general = (length(find(psth_general.bootsTrapTest == -1)) / length(psth_general.bootsTrapTest))*100;

% Visual Object
psth_visualObject = spikesPsth(pulses.timestampsOn{3}(indexes == 1)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_visualObject.responsecurveSmooth,1)
    psth_visualObject.responseZ(ii,:) = psth_visualObject.responsecurveZSmooth(ii,:) - ...
        mean(psth_visualObject.responsecurveSmooth(ii,win_Z))./std(psth_visualObject.responsecurveSmooth(ii,win_Z));
end
psth_visualObject.peakResponse = nanmean(psth_visualObject.responsecurve(:,win),2);
psth_visualObject.peakResponseZ = nanmean(psth_visualObject.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_visualObject.responsecurveZSmooth(find(~isnan(psth_visualObject.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_visualObject.peakResponseZ(find(~isnan(psth_visualObject.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_visualObject.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_visualObject = (length(find(psth_visualObject.bootsTrapTest == 1)) / length(psth_visualObject.bootsTrapTest))*100;
negative_responsive_cells_visualObject = (length(find(psth_visualObject.bootsTrapTest == -1)) / length(psth_visualObject.bootsTrapTest))*100;

% Spatial
psth_spatial = spikesPsth(pulses.timestampsOn{3}(indexes == 2)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_spatial.responsecurveSmooth,1)
    psth_spatial.responseZ(ii,:) = psth_spatial.responsecurveZSmooth(ii,:) - ...
        mean(psth_spatial.responsecurveSmooth(ii,win_Z))./std(psth_spatial.responsecurveSmooth(ii,win_Z));
end
psth_spatial.peakResponse = nanmean(psth_spatial.responsecurve(:,win),2);
psth_spatial.peakResponseZ = nanmean(psth_spatial.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_spatial.responsecurveZSmooth(find(~isnan(psth_spatial.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_spatial.peakResponseZ(find(~isnan(psth_spatial.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_spatial.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_spatial = (length(find(psth_spatial.bootsTrapTest == 1)) / length(psth_spatial.bootsTrapTest))*100;
negative_responsive_cells_spatial = (length(find(psth_spatial.bootsTrapTest == -1)) / length(psth_spatial.bootsTrapTest))*100;

% Abstract
psth_abstract = spikesPsth(pulses.timestampsOn{3}(indexes == 3)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_abstract.responsecurveSmooth,1)
    psth_abstract.responseZ(ii,:) = psth_abstract.responsecurveZSmooth(ii,:) - ...
        mean(psth_abstract.responsecurveSmooth(ii,win_Z))./std(psth_abstract.responsecurveSmooth(ii,win_Z));
end
psth_abstract.peakResponse = nanmean(psth_abstract.responsecurve(:,win),2);
psth_abstract.peakResponseZ = nanmean(psth_abstract.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_abstract.responsecurveZSmooth(find(~isnan(psth_abstract.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_abstract.peakResponseZ(find(~isnan(psth_abstract.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_abstract.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_abstract = (length(find(psth_abstract.bootsTrapTest == 1)) / length(psth_abstract.bootsTrapTest))*100;
negative_responsive_cells_abstract = (length(find(psth_abstract.bootsTrapTest == -1)) / length(psth_abstract.bootsTrapTest))*100;

% Actions
psth_actions = spikesPsth(pulses.timestampsOn{3}(indexes == 4)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_actions.responsecurveSmooth,1)
    psth_actions.responseZ(ii,:) = psth_actions.responsecurveZSmooth(ii,:) - ...
        mean(psth_actions.responsecurveSmooth(ii,win_Z))./std(psth_actions.responsecurveSmooth(ii,win_Z));
end
psth_actions.peakResponse = nanmean(psth_actions.responsecurve(:,win),2);
psth_actions.peakResponseZ = nanmean(psth_actions.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_actions.responsecurveZSmooth(find(~isnan(psth_actions.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_actions.peakResponseZ(find(~isnan(psth_actions.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_actions.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_action = (length(find(psth_actions.bootsTrapTest == 1)) / length(psth_actions.bootsTrapTest))*100;
negative_responsive_cells_action = (length(find(psth_actions.bootsTrapTest == -1)) / length(psth_actions.bootsTrapTest))*100;

% Familiar People
psth_familiar = spikesPsth(pulses.timestampsOn{3}(indexes == 5)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_familiar.responsecurveSmooth,1)
    psth_familiar.responseZ(ii,:) = psth_familiar.responsecurveZSmooth(ii,:) - ...
        mean(psth_familiar.responsecurveSmooth(ii,win_Z))./std(psth_familiar.responsecurveSmooth(ii,win_Z));
end
psth_familiar.peakResponse = nanmean(psth_familiar.responsecurve(:,win),2);
psth_familiar.peakResponseZ = nanmean(psth_familiar.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_familiar.responsecurveZSmooth(find(~isnan(psth_familiar.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_familiar.peakResponseZ(find(~isnan(psth_familiar.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_familiar.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_familiar = (length(find(psth_familiar.bootsTrapTest == 1)) / length(psth_familiar.bootsTrapTest))*100;
negative_responsive_cells_familiar = (length(find(psth_familiar.bootsTrapTest == -1)) / length(psth_familiar.bootsTrapTest))*100;

%% ======= ARRIBA vs ABAJO =============

arriba = strcmpi(words,'Arriba');
arriba_index = find(arriba(:,1) == 1);

psth_arriba = spikesPsth(pulses.timestampsOn{3}(arriba_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);


abajo = strcmpi(words,'Abajo');
abajo_index = find(abajo(:,1) == 1);

psth_abajo = spikesPsth(pulses.timestampsOn{3}(abajo_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);



%% ======= IZQUIERDA vs DERECHA =============

izquierda = strcmpi(words,'Izquierda');
izquierda_index = find(izquierda(:,1) == 1);

psth_izquierda = spikesPsth(pulses.timestampsOn{3}(izquierda_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);


derecha = strcmpi(words,'Derecha');
derecha_index = find(derecha(:,1) == 1);

psth_derecha = spikesPsth(pulses.timestampsOn{3}(derecha_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);




%% ========== REACTION TIME (-4 seconds) ===================

% General
psth_general = spikesPsth([pulses.timestampsOn{3}'-4],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

win_resp = [0.8 1.1];
ts_psth = psth_general.timestamps;
win_Z = find(ts_psth >= -1 & ts_psth <= -0.1);
win = find(ts_psth >= win_resp(1) & ts_psth <= win_resp(2));

for ii = 1:size(psth_general.responsecurveSmooth,1)
    psth_general.responseZ(ii,:) = psth_general.responsecurveZSmooth(ii,:) - ...
        mean(psth_general.responsecurveSmooth(ii,win_Z))./std(psth_general.responsecurveSmooth(ii,win_Z));
end
psth_general.peakResponse = nanmean(psth_general.responsecurve(:,win),2);
psth_general.peakResponseZ = nanmean(psth_general.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_general.responsecurveZSmooth(find(~isnan(psth_general.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_general.peakResponseZ(find(~isnan(psth_general.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_general.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_general = (length(find(psth_general.bootsTrapTest == 1)) / length(psth_general.bootsTrapTest))*100;
negative_responsive_cells_general = (length(find(psth_general.bootsTrapTest == -1)) / length(psth_general.bootsTrapTest))*100;

% Visual Object
psth_visualObject = spikesPsth(pulses.timestampsOn{3}(indexes == 1)'-4,'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_visualObject.responsecurveSmooth,1)
    psth_visualObject.responseZ(ii,:) = psth_visualObject.responsecurveZSmooth(ii,:) - ...
        mean(psth_visualObject.responsecurveSmooth(ii,win_Z))./std(psth_visualObject.responsecurveSmooth(ii,win_Z));
end
psth_visualObject.peakResponse = nanmean(psth_visualObject.responsecurve(:,win),2);
psth_visualObject.peakResponseZ = nanmean(psth_visualObject.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_visualObject.responsecurveZSmooth(find(~isnan(psth_visualObject.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_visualObject.peakResponseZ(find(~isnan(psth_visualObject.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_visualObject.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_visualObject = (length(find(psth_visualObject.bootsTrapTest == 1)) / length(psth_visualObject.bootsTrapTest))*100;
negative_responsive_cells_visualObject = (length(find(psth_visualObject.bootsTrapTest == -1)) / length(psth_visualObject.bootsTrapTest))*100;

% Spatial
psth_spatial = spikesPsth(pulses.timestampsOn{3}(indexes == 2)'-4,'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_spatial.responsecurveSmooth,1)
    psth_spatial.responseZ(ii,:) = psth_spatial.responsecurveZSmooth(ii,:) - ...
        mean(psth_spatial.responsecurveSmooth(ii,win_Z))./std(psth_spatial.responsecurveSmooth(ii,win_Z));
end
psth_spatial.peakResponse = nanmean(psth_spatial.responsecurve(:,win),2);
psth_spatial.peakResponseZ = nanmean(psth_spatial.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_spatial.responsecurveZSmooth(find(~isnan(psth_spatial.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_spatial.peakResponseZ(find(~isnan(psth_spatial.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_spatial.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_spatial = (length(find(psth_spatial.bootsTrapTest == 1)) / length(psth_spatial.bootsTrapTest))*100;
negative_responsive_cells_spatial = (length(find(psth_spatial.bootsTrapTest == -1)) / length(psth_spatial.bootsTrapTest))*100;

% Abstract
psth_abstract = spikesPsth(pulses.timestampsOn{3}(indexes == 3)'-4,'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_abstract.responsecurveSmooth,1)
    psth_abstract.responseZ(ii,:) = psth_abstract.responsecurveZSmooth(ii,:) - ...
        mean(psth_abstract.responsecurveSmooth(ii,win_Z))./std(psth_abstract.responsecurveSmooth(ii,win_Z));
end
psth_abstract.peakResponse = nanmean(psth_abstract.responsecurve(:,win),2);
psth_abstract.peakResponseZ = nanmean(psth_abstract.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_abstract.responsecurveZSmooth(find(~isnan(psth_abstract.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_abstract.peakResponseZ(find(~isnan(psth_abstract.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_abstract.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_abstract = (length(find(psth_abstract.bootsTrapTest == 1)) / length(psth_abstract.bootsTrapTest))*100;
negative_responsive_cells_abstract = (length(find(psth_abstract.bootsTrapTest == -1)) / length(psth_abstract.bootsTrapTest))*100;

% Actions
psth_action = spikesPsth(pulses.timestampsOn{3}(indexes == 4)'-4,'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

for ii = 1:size(psth_action.responsecurveSmooth,1)
    psth_action.responseZ(ii,:) = psth_action.responsecurveZSmooth(ii,:) - ...
        mean(psth_action.responsecurveSmooth(ii,win_Z))./std(psth_action.responsecurveSmooth(ii,win_Z));
end
psth_action.peakResponse = nanmean(psth_action.responsecurve(:,win),2);
psth_action.peakResponseZ = nanmean(psth_action.responseZ(:,win),2);

figure;
imagesc_ranked(ts_psth,[],psth_action.responsecurveZSmooth(find(~isnan(psth_action.responsecurveZSmooth(:,1))),:),[-3 3],...
    psth_action.peakResponseZ(find(~isnan(psth_action.responsecurveZSmooth(:,1)))));
colormap jet;
hold on;
zmean = nanmean(zscore(psth_action.responsecurveZSmooth,[],2));
zmean = zmean - min(zmean); 
zmean = zmean/max(zmean) * (111-1) * std(zmean);
plot(ts_psth,zmean+10,'k','LineWidth',1);

positive_responsive_cells_action = (length(find(psth_action.bootsTrapTest == 1)) / length(psth_action.bootsTrapTest))*100;
negative_responsive_cells_action = (length(find(psth_action.bootsTrapTest == -1)) / length(psth_action.bootsTrapTest))*100;


%% ============== REACTION TIMES ======================

reaction_times = cell2mat(words(:,3));

reaction_times_visualObject = reaction_times(indexes == 1);
reaction_times_spatial = reaction_times(indexes == 2);
reaction_times_abstract = reaction_times(indexes == 3);
reaction_times_action = reaction_times(indexes == 4);







%% FLOR 

flor = strcmpi(words,'Flor');
flor_index = find(flor(:,1) == 1);

psth_flor = spikesPsth(pulses.timestampsOn{3}(flor_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);

%% IZQUIERDA
izquierda = strcmpi(words,'Izquierda');
izquierda_index = find(izquierda(:,1) == 1);

psth_izquierda = spikesPsth(pulses.timestampsOn{3}(izquierda_index)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',0,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2], 'binSize', 0.01, 'win_Z',[-3 1]);



%% ========= ANALOG INPUT 2 ==============
if ~isempty(dir('semantic_words_2.mat'))
    targetFile = dir('semantic_words_2.mat');
    load(targetFile.name);
end

% General
psth_general_1 = spikesPsth(pulses.timestampsOn{2}','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 2], 'binSize', 0.01, 'win_Z',[-3 1]);












end