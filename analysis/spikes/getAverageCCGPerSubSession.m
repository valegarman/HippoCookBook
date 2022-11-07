function [averageCCGSubSessions] = getAverageCCGPerSubSession(varargin)
% [averageCCGSubSession] = getAverageCCGPerSubSession(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% for each subfolder.
%
% <OPTIONALS>
% spikes        buzcode spikes structure, if not provided tries loadSpikes.
% basepath      By default pwd.
% binSize       In seconds, default, 0.005.
% winSize       In seconds, default, 0.6.
% winSizePlot   Interval in seconds, default [-.3 .3] 
% plotOpt       Default true.
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
%
% OUTPUTS
% averageCCG
%
% Pablo Abad 2022 based on Manu-BuzsakiLab 2021

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
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'includeIntervals',[0 Inf],@isnumeric);
addParameter(p,'winIndex',[-.01 .01],@isnumeric);
addParameter(p,'interp0',[-.01 .01],@isnumeric);
addParameter(p,'useBrainRegions',true,@islogical);
addParameter(p,'useDistinctShanks',true,@islogical);

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
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;
includeIntervals = p.Results.includeIntervals;
winIndex = p.Results.winIndex;
interp0 = p.Results.interp0;
useBrainRegions = p.Results.useBrainRegions;
useDistinctShanks = p.Results.useDistinctShanks;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.averageCCGSubSessions.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Average CCG already computed! Loading file...');
    load(targetFile.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

try
    targetFile = dir('*.MergePoints.events.mat');
    if ~isempty(targetFile)
        disp('MergePoints detected. Loading file...');
        load(targetFile.name);
    end
catch
    warning('MergePoints not found. Skipping ...');
    return;
end

if skipStimulationPeriods
    try
        try 
            targetFile = dir('*optogeneticPulses*');
            optogenetic_responses = importdata(targetFile.name);
        catch
            warning('Could not open optogeneticPulses file... trying to open optogeneticResponses...');
            optogenetic_responses = getOptogeneticResponse;
        end
        excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
    catch
        warning('Skip stimulation periods not possible...');
    end
end

if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end

averageCCGSubSessions = [];

% Now we call getAverageCCG with includeIntervals

for ii = 1:length(MergePoints.foldernames)
    includeIntervals = MergePoints.timestamps(ii,:);
    averageCCG = getAverageCCG('includeIntervals',includeIntervals,'force',true,'saveMat',false,'saveFig',false);
    
    averageCCGSubSessions.(MergePoints.foldernames{ii}) = averageCCG;
end

% Save output
if saveMat
    disp('Saving results...');
    session = loadSession;
    filename = session.general.name;
    save([filename '.averageCCGSubSessions.cellinfo.mat'],'averageCCGSubSessions');
end


% Plotting general figure
if plotOpt
    figure;
    set(gcf,'Position',[200 -500 2500 1200]);
    brColors = hsv(length(fields(averageCCGSubSessions)));
    fldnames = cell(length(fields(averageCCGSubSessions))*2,1);
    fldnames(2:2:end) = fields(averageCCGSubSessions);
    fldnames(1:2:end) = {' '};
    for ii = 1:length(fields(averageCCGSubSessions))
        t_ccg = averageCCGSubSessions.(MergePoints.foldernames{ii}).timestamps;
        win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
        plotFill(t_ccg,averageCCGSubSessions.(MergePoints.foldernames{ii}).ZmedianCCG,'faceAlpha',0.1,'color',brColors(ii,:));
        hold on;

    end
    legend(fldnames)
    
    if saveFig
        saveas(gcf,['SummaryFigures\averageCCGSubSessions.png']);
    end
end

close all;
cd(prevPath);
end