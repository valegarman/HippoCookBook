function [spikeCCGchange] = getSpikeCCGchange(referenceNeuron,intervals,varargin)
% [spikeEpochsCCG] = getSpikeCCGchange(varargin)
%
% Computes CCG between a neuron and the remaining cells in different epochs
% and change statistics
%
% referenceNeuron       UID of the reference neuron. If empty, gets neuron
%                           from spikeTriggeredPulse structure
% intervals             Intervals to compute the CCG in. If empty, gets pre
%                           and post stimulation intervals in the
%                           spikeTriggeredPulse structure
%
% <OPTIONALS>
% spikes                buzcode spikes structure, if not provided tries loadSpikes.
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
% iluminated
% OUTPUTS
% averageCCG
%
% Manu-BuzsakiLab 2021

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
addParameter(p,'winIndex',[-.01 .01],@isnumeric);
addParameter(p,'interp0',[-.01 .01],@isnumeric);
addParameter(p,'useBrainRegions',true,@islogical);
addParameter(p,'useDistinctShanks',true,@islogical);
% addParameter(p,'useCellType',true,@islogical); work in progress
addParameter(p,'restrict_ints',[0 Inf],@isnumeric);
addParameter(p,'save_as','spikeCCGchange',@ischar);
addParameter(p,'zwin',[-0.3 -0.15],@isnumeric);

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
winIndex = p.Results.winIndex;
interp0 = p.Results.interp0;
useBrainRegions = p.Results.useBrainRegions;
useDistinctShanks = p.Results.useDistinctShanks;
% useCellType = p.Results.useCellType;
restrict_ints = p.Results.restrict_ints;
save_as = p.Results.save_as;
zwin = p.Results.zwin;

% Deal with inputs
prevPath = pwd;
cd(basepath);

session = loadSession;

targetFile = dir('*.spikeCCGchange.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Spikes epoch CCG already computed! Loading file...');
    load(targetFile.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if any(restrict_ints ~= [0 Inf])
    warning('Restricting analysis for intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},restrict_ints);
        spikes.times{ii} = spikes.times{ii}(status);
    end 
end

% do ccg
intervalCCG = [];
for ii = 1:size(intervals,1)
    spikesTimesIntervals = spikes.times;
    for jj = 1:length(spikesTimesIntervals)
        spikesTimesIntervals{jj}(~InIntervals(spikesTimesIntervals{jj},intervals(ii,:))) = [];
    end

    [allCcg, t_ccg] = CCG(spikesTimesIntervals,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);

    intervalCCG_temp = squeeze(allCcg(:,referenceNeuron,:));
    intervalCCG_temp(:,referenceNeuron) = NaN;

    intervalCCG = cat(3,intervalCCG,intervalCCG_temp);
    intervalCCGZ(:,:,ii) = zscore_win(squeeze(intervalCCG(:,:,ii)),InIntervals(t_ccg,zwin)); 
end

% keyboard;
% shank3 = spikes.shankID==3;
% shank1 = spikes.shankID==1;
% drivenCells_shank3 = spikes.shankID==3 & uLEDResponses.drivenCells';
% notDrivenCells_shank3 = spikes.shankID==3 & ~uLEDResponses.drivenCells';

cofiringIndex = squeeze(median(intervalCCG(InIntervals(t_ccg, winIndex),:,:),1));
cofiringIndexZ = squeeze(median(intervalCCGZ(InIntervals(t_ccg, winIndex),:,:),1));

spikeCCGchange.cofiringIndex = cofiringIndex;
spikeCCGchange.cofiringIndexZ = cofiringIndexZ;
spikeCCGchange.intervalCCG = intervalCCG;
spikeCCGchange.intervalCCGZ = intervalCCGZ;
spikeCCGchange.winIndex = winIndex;
spikeCCGchange.timestamps = t_ccg;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename, '.', save_as, '.cellinfo.mat'],'spikeCCGchange');
end

close all;
cd(prevPath);
end