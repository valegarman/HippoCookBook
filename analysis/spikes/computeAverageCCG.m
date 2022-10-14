function [averageCCG] = computeAverageCCG(varargin)
% [averageCCG] = getAverageCCG(varargin)
%
% Computes Psth and a several statistical measures of the cell responses,
% taking into account the brain regions and comparing between all of them.
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
% Pablo Abad 2022, based on getAverageCCG by Manuel Valero 2022.

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
addParameter(p,'useBrainRegions',true,@islogical);
addParameter(p,'useDistinctShanks',true,@islogical);

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
useBrainRegions = p.Results.useBrainRegions;
useDistinctShanks = p.Results.useDistinctShanks;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.averageRegionCCG.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Average CCG per Region already computed! Loading file...');
    load(targetFile.name);
    return
end

% Possible layers per region (and possible regions)
CA1 = {'CA1so','CA1sp','CA1sr','CA1slm'};
CA2 = {'CA2so','CA2sp','CA2sr','CA2slm'};
CA3 = {'CA3so','CA3sp','CA3sr','CA3slm','CA3slu'};

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
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

%% Compute CCG
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

%% Brain Regions
brainRegionCCG = [];
if useBrainRegions && exist([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat'])
    disp('Computing CCG by brain region');
    session = loadSession;
    load([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat']);
    if isfield(session,'brainRegions')
        isCA1 = find(ismember(cell_metrics.brainRegion,CA1));
        isCA2 = find(ismember(cell_metrics.brainRegion,CA2));
        isCA3 = find(ismember(cell_metrics.brainRegion,CA3));
        for ii = 1:length(isCA1)
            cell_metrics.brainRegion{isCA1(ii)} = 'CA1';
        end
        for ii = 1:length(isCA2)
            cell_metrics.brainRegion{isCA2(ii)} = 'CA2';
        end
        for ii = 1:length(isCA3)
            cell_metrics.brainRegion{isCA3(ii)} = 'CA3';
        end

        efields = unique(cell_metrics.brainRegion);

        for ii = 1:length(efields)
            % Cells in same region
            cellsInRegion = ismember(cell_metrics.brainRegion,efields{ii});
            for jj = 1:length(spikes.times)
                cellsID =  indCell(indCell ~= jj & cellsInRegion);

                ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
                ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG

                ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
                ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG
            end

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

            brainRegionCCG.([efields{ii} '_medianCCG']) = ccMedian;
            brainRegionCCG.([efields{ii} '_ZmedianCCG']) = ccZMedian;
            brainRegionCCG.([efields{ii} '_meanCCG']) = ccMean;
            brainRegionCCG.([efields{ii} '_ZmeanCCG']) = ccZMean;
            brainRegionCCG.([efields{ii} '_ccgIndex']) = ccgIndex;
            brainRegionCCG.binSize = binSize;
            brainRegionCCG.winSize = winSize;
            brainRegionCCG.timestamps = t_ccg;
            brainRegionCCG.excludeIntervals = excludeIntervals;
            brainRegionCCG.winIndex = winIndex;
            brainRegionsCcgIndex.(efields{ii}) = ccgIndex;
        end
    else
        warning('Brain regions have not been defined yet...');
    end
     % CCGIndex per region
    efields = fieldnames(brainRegionsCcgIndex);
    for jj = 1:size(spikes.UID,2)
        for ii = 1:length(efields)
            ccgIndexRegion(jj,ii) = brainRegionsCcgIndex.(efields{ii})(jj);
        end
    end
    brainRegionCCG.ccgIndexRegion = ccgIndexRegion;
    brainRegionCCG.absCcgIndexRegion = abs(ccgIndexRegion);
    brainRegionCCG.listOfRegions = efields;
    brainRegionCCG.listOfRegionsID = 1:length(efields);
    brainRegionCCG.brainRegionsCcgIndex = brainRegionsCcgIndex;
    
    averageCCG.brainRegionCCG = brainRegionCCG;
end


if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.averageRegionCCG.cellinfo.mat'],'averageCCG');
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
        title([num2str(jj) ,'-',cell_metrics.brainRegion{jj}],'FontWeight','normal','FontSize',10);

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
    
    % by BrainRegion
    if isfield(averageCCG,'brainRegionCCG')
       
        % brainRegions colors
        brColors = hsv(length(averageCCG.brainRegionCCG.listOfRegionsID));
        
        figure;
        set(gcf,'Position',[200 -500 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            subplot(7,ceil(size(spikes.UID,2)/7),jj);
            hold on
            b = bar(averageCCG.brainRegionCCG.listOfRegionsID, averageCCG.brainRegionCCG.absCcgIndexRegion(jj,:));
            b.FaceColor = 'flat';
            b.CData = brColors;
            
            set(gca,'TickDir','out','XTick',averageCCG.brainRegionCCG.listOfRegionsID,'XTickLabel',averageCCG.brainRegionCCG.listOfRegions,'XTickLabelRotation',45);
            
            title([num2str(jj),'-',cell_metrics.brainRegion{jj}],'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('CCG Index');
            end
        end
        saveas(gcf,['SummaryFigures\CCGAvgPerRegion.png']); 
        
        % All cells per region
        regions = averageCCG.brainRegionCCG.listOfRegions;
        
        for ii = 1:length(regions)
            for jj = 1:length(regions)
                figure,
                set(gcf,'Position',[200 -500 2500 1200]);
                for kk = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),kk);
                    inRegion1 = ismember(cell_metrics.brainRegion,regions{ii});
                    inRegion2 = ismember(cell_metrics.brainRegion,regions{jj});
                    if inRegion1(kk) == 1
                        cc = zscore(squeeze(allCcg(:,kk,inRegion2))',[],2);
                        imagesc(t_ccg,1:length(inRegion2),cc)
                        hold on;
                        set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                        zmean = mean(zscore(squeeze(allCcg(:,kk,inRegion2))',[],2));
                        zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                        plot(t_ccg, zmean,'k','LineWidth',2);
                        title([num2str(kk) ,'-',cell_metrics.brainRegion{kk},'-',regions{jj}],'FontWeight','normal','FontSize',10);
                        if kk == 1
                            ylabel('Cell');
                        elseif kk == size(spikes.UID,2)
                            xlabel('Time (s)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                end
                saveas(gcf,['SummaryFigures\allCellsAverageCCG_',regions{ii},'_',regions{jj},'.png']);
            end
        end
        
        % Plot mean CCG per area and different area
        regions = averageCCG.brainRegionCCG.listOfRegions;
        figure;
        set(gcf,'Position',[200 -500 2500 1200]);
        counter = 0;
        for ii = 1:length(regions)
            for jj = 1:length(regions)
                counter = counter + 1;
                subplot(length(regions),2,counter)
                inRegion = ismember(cell_metrics.brainRegion,regions{ii});
                zmeanCCG = averageCCG.brainRegionCCG.([regions{jj},'_ZmeanCCG'])(inRegion,:);
                plotFill(t_ccg,zmeanCCG,'color',[0 0 0]);
                title(['Region: ' , regions{ii}, ' & ', regions{jj}])
                
                ylabel('ZmeanCCG');
                xlabel('Time (s)');
            end
        end
        saveas(gcf,['SummaryFigures\CCGRegionAverage.png']);
    end
end

cd(prevPath);
end