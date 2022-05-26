
function [averageCCG] = getAverageCCG(varargin)
% [averageCCG] = getAverageCCG(varargin)
%
% Computes Psth and a several statistical measures of the cell responses.
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
%
% OUTPUTS
% averageCCG
%
% Manu-BuzsakiLab 2021

% Parse options
p = inputParser;
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'binSize',0.005,@isnumeric);
addParameter(p,'winSize',0.6,@isnumeric);
addParameter(p,'plotOpt',true,@islogical);
addParameter(p,'winSizePlot',[-.3 .3],@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'winIndex',[-.01 .01],@isnumeric);
addParameter(p,'interp0',[-.01 .01],@isnumeric);

parse(p, varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
winSizePlot = p.Results.winSizePlot;
plotOpt = p.Results.plotOpt;
saveMat = p.Results.saveMat;
force = p.Results.force;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;
winIndex = p.Results.winIndex;
interp0 = p.Results.interp0;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.averageCCG.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Average CCG already computed! Loading file...');
    load(targetFile.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if skipStimulationPeriods
    try
        optogenetic_responses = getOptogeneticResponse;
    catch
        warning('Skip stimulation periods not possible...');
    end
end
excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end
% do ccg
[allCcg, t_ccg] = CCG(spikes.times,[],'binSize',binSize,'duration',winSize);
indCell = [1:size(allCcg,2)];
for jj = 1 : length(spikes.times)
    ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,indCell(indCell~=jj))),2); % zCCG
    ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2)); % zCCG
    
    ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,indCell(indCell~=jj))),2); % zCCG
    ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2)); % zCCG
end

% interpolate value in 0
if interp0
    artifactSamples = find(t_ccg == 0);
    x_axis = 1:length(t_ccg);
    x_axis(artifactSamples) = [];
    for jj = 1:size(ccMedian,1)
        ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
        ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
        ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
        ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
    end
end

win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
ccgIndex = median(ccZMedian(:,win),2);

averageCCG.medianCCG = ccMedian;
averageCCG.ZmedianCCG = ccZMedian;
averageCCG.meanCCG = ccMean;
averageCCG.ZmeanCCG = ccZMean;
averageCCG.binSize = binSize;
averageCCG.winSize = winSize;
averageCCG.timestamps = t_ccg;
averageCCG.excludeIntervals = excludeIntervals;
averageCCG.ccgIndex = ccgIndex;
averageCCG.winIndex = winIndex;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.averageCCG.cellinfo.mat'],'averageCCG');
end

if plotOpt
    % all cells
    figure;
    set(gcf,'Position',[200 -500 2500 1200]);
    for jj = 1:size(spikes.UID,2)
        % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj);
        cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
        imagesc(t_ccg,1:max(indCell)-1,cc)
        set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
        hold on
        zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
        zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
        plot(t_ccg, zmean,'k','LineWidth',2);
        xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
        title(num2str(jj),'FontWeight','normal','FontSize',10);

        if jj == 1
            ylabel('Cell');
        elseif jj == size(spikes.UID,2)
            xlabel('Time (s)');
        else
            set(gca,'YTick',[],'XTick',[]);
        end
    end
    saveas(gcf,['SummaryFigures\allCellsAverageCCG.png']); 
    
    % grand mean
    figure
    imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG.ZmeanCCG,1)],...
        averageCCG.ZmeanCCG); caxis([-3 3]); colormap(jet);
    set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
    title('Grand CCG average','FontWeight','normal','FontSize',10);
    saveas(gcf,['SummaryFigures\grandCCGAverage.png']);
end

cd(prevPath);
end