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
% Manu-BuzsakiLab 2021
% Manu-NCL 2024, improving Z-score, regions, cell-types, etc

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
addParameter(p,'interp0',false,@islogical);
addParameter(p,'useBrainRegions',true,@islogical);
addParameter(p,'useDistinctShanks',true,@islogical);
% addParameter(p,'useCellType',true,@islogical); work in progress
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','averageCCG',@ischar);
addParameter(p,'win_Z',[-0.3 -0.15],@isnumeric);
addParameter(p,'orderOfccZMedianMap',[]);

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
% useCellType = p.Results.useCellType;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;
win_Z = p.Results.win_Z;
orderOfccZMedianMap = p.Results.orderOfccZMedianMap;

% Deal with inputs
prevPath = pwd;
cd(basepath);

% Load session
session = loadSession();

targetFile = dir('*.averageCCG.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Average CCG already computed! Loading file...');
    load(targetFile.name);
    return
end

ints = [];
session = loadSession;
if isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = 'averageCCG_post';
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_baseline
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [0 session.epochs{ii}.startTime];
            warning('Epoch with manipulations found! Restricting analysis to baseline interval!');
        end
    end
    if isempty(ints)
        ints = [0 Inf];
    end
else
    ints = [0 Inf];
end
restrict_ints = IntersectIntervals([ints; restrict_to]);

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

if includeIntervals(1)>0 || includeIntervals(2)<Inf
    warning('Including intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},includeIntervals);
        spikes.times{ii} = spikes.times{ii}(status);
    end
end

if any(restrict_ints ~= [0 Inf])
    warning('Restricting analysis for intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},restrict_ints);
        spikes.times{ii} = spikes.times{ii}(status);
    end 
end

% do ccg
[allCcg, t_ccg] = CCG(spikes.times,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);
averageCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, ones(1,size(allCcg,2)), orderOfccZMedianMap);
averageCCG.excludeIntervals = excludeIntervals;

brainRegionCCG = [];
if useBrainRegions && exist([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat'])
    disp('Computing CCG by brain region');
    cell_metrics = loadCellMetrics;
    session = loadSession;
    
    if isfield(session,'brainRegions')
        efields = fieldnames(session.brainRegions);
        for ii =  1: length(efields)
            cellsInRegion = ismember(cell_metrics.brainRegion,efields{ii});
            avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0,cellsInRegion);
            
            brainRegionCCG.([efields{ii} '_medianCCG']) = avgCCG.medianCCG;
            brainRegionCCG.([efields{ii} '_ZmedianCCG']) = avgCCG.ZmedianCCG;
            brainRegionCCG.([efields{ii} '_meanCCG']) = avgCCG.meanCCG;
            brainRegionCCG.([efields{ii} '_ZmeanCCG']) = avgCCG.ZmeanCCG;
            brainRegionCCG.([efields{ii} '_ccgIndex']) = avgCCG.ccgIndex;
            brainRegionCCG.([efields{ii} '_ccgMeanIndex']) = avgCCG.ccgMeanIndex;
            brainRegionCCG.([efields{ii} '_ccgAbsMeanIndex']) = avgCCG.ccgAbsMeanIndex;
            brainRegionCCG.([efields{ii} '_ccgMedianIndex']) = avgCCG.ccgMedianIndex;
            brainRegionCCG.([efields{ii} '_ccgAbsMedianIndex']) = avgCCG.ccgAbsMedianIndex;
            brainRegionCCG.([efields{ii} '_ccZMedianMap']) = avgCCG.ccZMedianMap;

            brainRegionCCG.binSize = avgCCG.binSize;
            brainRegionCCG.winSize = avgCCG.winSize;
            brainRegionCCG.timestamps = avgCCG.timestamps;
            brainRegionCCG.winIndex = avgCCG.winIndex;
            brainRegionCCG.win_Z = avgCCG.win_Z;

            brainRegionCCG.excludeIntervals = excludeIntervals;
            brainRegionsCcgIndex.(efields{ii}) = avgCCG.ccgIndex;
        end
        
        % for CA1
        cellsInRegion = ismember(cell_metrics.brainRegion,'CA1') | ismember(cell_metrics.brainRegion,'CA1sp')...
            | ismember(cell_metrics.brainRegion,'CA1so') | ismember(cell_metrics.brainRegion,'CA1slm') | ismember(cell_metrics.brainRegion,'CA1sr');
        avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, cellsInRegion);
        

        brainRegionCCG.('CA1_medianCCG') = avgCCG.medianCCG;
        brainRegionCCG.('CA1_ZmedianCCG') = avgCCG.ZmedianCCG;
        brainRegionCCG.('CA1_meanCCG') = avgCCG.meanCCG;
        brainRegionCCG.('CA1_ZmeanCCG') = avgCCG.ZmeanCCG;
        brainRegionCCG.('CA1_ccgIndex') = avgCCG.ccgIndex;
        brainRegionsCcgIndex.('CA1') = avgCCG.ccgIndex;
        brainRegionCCG.('CA1_ccgMeanIndex') = avgCCG.ccgMeanIndex;
        brainRegionCCG.('CA1_ccgAbsMeanIndex') = avgCCG.ccgAbsMeanIndex;
        brainRegionCCG.('CA1_ccgMedianIndex') = avgCCG.ccgMedianIndex;
        brainRegionCCG.('CA1_ccgAbsMedianIndex') = avgCCG.ccgAbsMedianIndex;
        brainRegionCCG.('CA1_ccZMedianMap') = avgCCG.ccZMedianMap;
        
        % for HPC
        cellsInRegion = ismember(cell_metrics.brainRegion,'CA1') | ismember(cell_metrics.brainRegion,'CA1sp')...
            | ismember(cell_metrics.brainRegion,'CA1so') | ismember(cell_metrics.brainRegion,'CA1slm') | ismember(cell_metrics.brainRegion,'CA1sr') ...
            | ismember(cell_metrics.brainRegion,'CA3') | ismember(cell_metrics.brainRegion,'CA3slm') | ismember(cell_metrics.brainRegion,'CA3slu') ...
            | ismember(cell_metrics.brainRegion,'CA3so') | ismember(cell_metrics.brainRegion,'CA3sp') | ismember(cell_metrics.brainRegion,'CA3sr') ...
            | ismember(cell_metrics.brainRegion,'CA2') | ismember(cell_metrics.brainRegion,'CA2slm') | ismember(cell_metrics.brainRegion,'CA2so') ...
            | ismember(cell_metrics.brainRegion,'CA2sp') | ismember(cell_metrics.brainRegion,'CA2sr') | ismember(cell_metrics.brainRegion,'DG') ...
            | ismember(cell_metrics.brainRegion,'HIP') | ismember(cell_metrics.brainRegion,'HPF');
        
        avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, cellsInRegion);
        
        brainRegionCCG.('HIP_medianCCG') = avgCCG.medianCCG;
        brainRegionCCG.('HIP_ZmedianCCG') = avgCCG.ZmedianCCG;
        brainRegionCCG.('HIP_meanCCG') = avgCCG.meanCCG;
        brainRegionCCG.('HIP_ZmeanCCG') = avgCCG.ZmeanCCG;
        brainRegionCCG.('HIP_ccgIndex') = avgCCG.ccgIndex;
        brainRegionsCcgIndex.('HIP') = avgCCG.ccgIndex;
        brainRegionCCG.('HIP_ccgMeanIndex') = avgCCG.ccgMeanIndex;
        brainRegionCCG.('HIP_ccgAbsMeanIndex') = avgCCG.ccgAbsMeanIndex;
        brainRegionCCG.('HIP_ccgMedianIndex') = avgCCG.ccgMedianIndex;
        brainRegionCCG.('HIP_ccgAbsMedianIndex') = avgCCG.ccgAbsMedianIndex;
        brainRegionCCG.('HIP_ccZMedianMap') = avgCCG.ccZMedianMap;

        % for CA1 pyramidal cells
        cellsInRegion = ismember(cell_metrics.brainRegion,'CA1') | ismember(cell_metrics.brainRegion,'CA1sp')...
            | ismember(cell_metrics.brainRegion,'CA1so') | ismember(cell_metrics.brainRegion,'CA1slm') | ismember(cell_metrics.brainRegion,'CA1sr');
        isPyramidal = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
        avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, cellsInRegion & isPyramidal);

        brainRegionCCG.('CA1pyr_medianCCG') = avgCCG.medianCCG;
        brainRegionCCG.('CA1pyr_ZmedianCCG') = avgCCG.ZmedianCCG;
        brainRegionCCG.('CA1pyr_meanCCG') = avgCCG.meanCCG;
        brainRegionCCG.('CA1pyr_ZmeanCCG') = avgCCG.ZmeanCCG;
        brainRegionCCG.('CA1pyr_ccgIndex') = avgCCG.ccgIndex;
        brainRegionsCcgIndex.('CA1pyr') = avgCCG.ccgIndex;
        brainRegionCCG.('CA1pyr_ccgMeanIndex') = avgCCG.ccgMeanIndex;
        brainRegionCCG.('CA1pyr_ccgAbsMeanIndex') = avgCCG.ccgAbsMeanIndex;
        brainRegionCCG.('CA1pyr_ccgMedianIndex') = avgCCG.ccgMedianIndex;
        brainRegionCCG.('CA1pyr_ccgAbsMedianIndex') = avgCCG.ccgAbsMedianIndex;
        brainRegionCCG.('CA1pyr_ccZMedianMap') = avgCCG.ccZMedianMap;
        
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
    else
        warning('Brain regions have not been defined yet...');
    end
end

session = loadSession;
shanksCCG = [];
if useDistinctShanks && length(session.extracellular.electrodeGroups.channels)>1 % if more than 1 shanks
    for ii =  1: length(session.extracellular.electrodeGroups.channels)
        cellsInShank = ismember(spikes.shankID,ii);
        avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, cellsInShank);
        

        shanksCCG.(['shank_' num2str(ii) '_medianCCG']) = avgCCG.medianCCG;
        shanksCCG.(['shank_' num2str(ii)  '_ZmedianCCG']) = avgCCG.ZmedianCCG;
        shanksCCG.(['shank_' num2str(ii)  '_meanCCG']) = avgCCG.meanCCG;
        shanksCCG.(['shank_' num2str(ii)  '_ZmeanCCG']) = avgCCG.ZmeanCCG;
        shanksCCG.(['shank_' num2str(ii)  '_ccgIndex']) = avgCCG.ccgIndex;
        shanksCCGIndex.(['shank_' num2str(ii)]) = avgCCG.ccgIndex;
        shanksCCG.(['shank_' num2str(ii)  '_ccgMeanIndex']) = avgCCG.ccgMeanIndex;
        shanksCCG.(['shank_' num2str(ii)  '_ccgAbsMeanIndex']) = avgCCG.ccgAbsMeanIndex;
        shanksCCG.(['shank_' num2str(ii)  '_ccgMedianIndex']) = avgCCG.ccgMedianIndex;
        shanksCCG.(['shank_' num2str(ii)  '_ccgAbsMedianIndex']) = avgCCG.ccgAbsMedianIndex;
        shanksCCG.(['shank_' num2str(ii)  '_ccZMedianMap']) = avgCCG.ccZMedianMap;
        shanksCCG.binSize = binSize;
        shanksCCG.winSize = winSize;
        shanksCCG.timestamps = t_ccg;
        shanksCCG.excludeIntervals = excludeIntervals;
        shanksCCG.winIndex = winIndex;
    end
    
    % CCGIndex per region
    clear ccgIndexShanks
    efields = fieldnames(shanksCCGIndex);
    for jj = 1:size(spikes.UID,2)
        for ii = 1:length(efields)
            ccgIndexShanks(jj,ii) = shanksCCGIndex.(efields{ii})(jj);
        end
    end
    shanksCCG.ccgIndexShanks = ccgIndexShanks;
    shanksCCG.absCcgIndexShanks = abs(ccgIndexShanks);
    shanksCCG.listOfShanks = efields;
    shanksCCG.listOfShankssID = 1:length(efields);
    shanksCCG.CcgIndexPerShank = shanksCCGIndex;

    averageCCG.shanksCCG = shanksCCG;
end

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename, '.', save_as, '.cellinfo.mat'],'averageCCG');
end

if plotOpt
    % all cells
    indCell = 1:size(allCcg,2);
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
    if saveFig
        saveas(gcf,['SummaryFigures\allCellsAverageCCG.png']); 
    end
    
    % grand mean
    figure
    imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG.ZmeanCCG,1)],...
        averageCCG.ZmeanCCG); caxis([-3 3]); colormap(jet);
    set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
    title('Grand CCG average','FontWeight','normal','FontSize',10);
    if saveFig
        saveas(gcf,['SummaryFigures\' save_as '.png']);
    end
    
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
            
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('CCG Index');
            end
        end
        if saveFig
            saveas(gcf,['SummaryFigures\', save_as,'PerRegion.png']); 
        end
    end
    
    % by shank
    if isfield(averageCCG,'shanksCCG')
        
        % brainRegions colors
        brColors = parula(length(averageCCG.shanksCCG.listOfShankssID));
        
        figure;
        set(gcf,'Position',[200 -500 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            subplot(7,ceil(size(spikes.UID,2)/7),jj);
            hold on
            b = bar(averageCCG.shanksCCG.listOfShankssID, averageCCG.shanksCCG.absCcgIndexShanks(jj,:));
            b.FaceColor = 'flat';
            b.CData = brColors;
            
            set(gca,'TickDir','out','XTick',averageCCG.shanksCCG.listOfShankssID,'XTickLabel',averageCCG.shanksCCG.listOfShanks,'XTickLabelRotation',45);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            
            plot(spikes.shankID(jj), averageCCG.shanksCCG.absCcgIndexShanks(jj,spikes.shankID(jj))+0.2,'.','MarkerSize',10,'Color',[.5 .5 .5]);
            text(spikes.shankID(jj), averageCCG.shanksCCG.absCcgIndexShanks(jj,spikes.shankID(jj))+0.4,'cell','Color',[.5 .5 .5]);
            
            if jj == 1
                ylabel('CCG Index');
            end
        end
        if saveFig
            saveas(gcf,['SummaryFigures\', save_as, 'PerShank.png']);
        end      
    end
    
end
% close all;
cd(prevPath);
end

function avgCCG = averageCCGheart(allCcg, t_ccg, win_Z, winIndex, interp0, include_cells, orderOfccZMedianMap)
    
    if nargin < 6
        include_cells = ones(1,size(allCcg,2));
    end

    if nargin < 7
        orderOfccZMedianMap = [];
    end

    if interp0
        artifactSamples = find(t_ccg == 0);
        x_axis = 1:length(t_ccg);
        x_axis(artifactSamples) = [];
    end
    
    status_winZ = InIntervals(t_ccg,win_Z);
    status_winIndex = InIntervals(t_ccg,winIndex);
    indCell = 1:size(allCcg,2);
    for jj = 1 : length(indCell)
        cellsID = indCell(indCell~=jj & include_cells);
        cc = squeeze(allCcg(:,jj,cellsID));
        % detecting of spike shadows, as negative peaks with 0 lag lower than
        % 1 std (the two is after accounting for the falling and rising step of
        % the shadow)
        center = ceil(size(cc,1)/2);
        spikes_shadow = diff(cc(center-1:center+1,:));
        spikes_shadow(2,:) = spikes_shadow(2,:) * -1;
        is_spikeshadow = sum(spikes_shadow) < -std(cc) * 2;
        cc(center,is_spikeshadow) = round(mean(cc([center-1 center+1], is_spikeshadow)));
        
        % stats
        ccMedian(jj,:) = median(cc,2); % 
        ccZMedian(jj,:) = median(zscore_win(cc,status_winZ),2);
        ccMean(jj,:) = mean(cc,2); % 
        ccZMean(jj,:) = mean(zscore_win(cc,status_winZ),2);
        if interp0
            ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
            ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
            ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
            ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
        end
        
        % index
        ccZ = zscore_win(cc,status_winZ,2);
        ccgIndexes(jj,:) = [mean(ccZMedian(jj,status_winIndex)) mean(ccZMean(jj,status_winIndex)) ...
            mean(median(abs(ccZ(status_winIndex,:)),2)) mean(mean(abs(ccZ(status_winIndex,:)),2))];

        % ccMap
        if ~isempty(orderOfccZMedianMap)
            orderOfMap2 = orderOfccZMedianMap(jj,:);
        else
            [~,orderOfMap2] =  sort(mean(ccZ(find(InIntervals(t_ccg,[-0.015 0.015])),:)));
        end
        
        if length(find(include_cells)) > 1
            ccZ = interpft(ccZ(:,orderOfMap2),100,2)';
            ccZMedianMap(jj,:,:) = ccZ;
            orderOfMap(jj,:) = nan(1, length(find(include_cells)));
            orderOfMap(jj,1:length(orderOfMap2)) = orderOfMap2;
        else
            ccZMedianMap(jj,:,:) = nan(size(t_ccg));
            orderOfMap(jj,:) = nan(1, length(find(include_cells)));
        end
        
    end
    
    avgCCG.medianCCG            = ccMedian;
    avgCCG.ZmedianCCG           = ccZMedian;
    avgCCG.meanCCG              = ccMean;
    avgCCG.ZmeanCCG             = ccZMean;
    avgCCG.ccgMedianIndex       = ccgIndexes(:,1);
    avgCCG.ccgMeanIndex         = ccgIndexes(:,2);
    avgCCG.ccgAbsMedianIndex    = ccgIndexes(:,4);
    avgCCG.ccgAbsMeanIndex      = ccgIndexes(:,4);
    avgCCG.ccgIndex             = ccgIndexes(:,1); % for retrocompatibility
    
    avgCCG.binSize              = median(diff(t_ccg));
    avgCCG.winSize              = max(t_ccg)*2;
    avgCCG.timestamps           = t_ccg;
    avgCCG.winIndex             = winIndex;
    avgCCG.win_Z                = win_Z;
    avgCCG.allCcg               = allCcg;
    avgCCG.ccZMedianMap         = ccZMedianMap;
    avgCCG.orderOfccZMedianMap  = orderOfMap;
end