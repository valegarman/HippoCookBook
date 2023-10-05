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
skipStimulationPeriods = p.Results.skipStimulationPeriods;
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
averageCCG.allCcg = allCcg;

brainRegionCCG = [];
if useBrainRegions && exist([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat'])
    disp('Computing CCG by brain region');
    session = loadSession;
    load([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat']);
    
    if isfield(session,'brainRegions')
        efields = fieldnames(session.brainRegions);
        for ii =  1: length(efields)
            cellsInRegion = ismember(cell_metrics.brainRegion,efields{ii});
            
            clear cellsID ccMedian ccZMedian ccMean ccZMean
            for jj = 1 : length(spikes.times)
                cellsID = indCell(indCell~=jj & cellsInRegion);
                if length(cellsID) == 1
                    ccMedian(jj,:) = squeeze(allCcg(:,jj,cellsID)); %
                    ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,cellsID))',[],2); % zCCG

                    ccMean(jj,:) = squeeze(allCcg(:,jj,cellsID)); % zCCG
                    ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,cellsID))',[],2); % zCCG
                else
                    ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
                    ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG

                    ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
                    ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG
                end
                
                ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
                ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG

                ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
                ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG
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
        
        % for CA1
        cellsInRegion = ismember(cell_metrics.brainRegion,'CA1') | ismember(cell_metrics.brainRegion,'CA1sp')...
            | ismember(cell_metrics.brainRegion,'CA1so') | ismember(cell_metrics.brainRegion,'CA1slm') | ismember(cell_metrics.brainRegion,'CA1sr');
        
        clear cellsID ccMedian ccZMedian ccMean ccZMean
        for jj = 1 : length(spikes.times)
            cellsID = indCell(indCell~=jj & cellsInRegion);

            ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
            ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG

            ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
            ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG
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

        brainRegionCCG.('CA1_medianCCG') = ccMedian;
        brainRegionCCG.('CA1_ZmedianCCG') = ccZMedian;
        brainRegionCCG.('CA1_meanCCG') = ccMean;
        brainRegionCCG.('CA1_ZmeanCCG') = ccZMean;
        brainRegionCCG.('CA1_ccgIndex') = ccgIndex;
        brainRegionsCcgIndex.('CA1') = ccgIndex;
        
        % for HPC
        cellsInRegion = ismember(cell_metrics.brainRegion,'CA1') | ismember(cell_metrics.brainRegion,'CA1sp')...
            | ismember(cell_metrics.brainRegion,'CA1so') | ismember(cell_metrics.brainRegion,'CA1slm') | ismember(cell_metrics.brainRegion,'CA1sr') ...
            | ismember(cell_metrics.brainRegion,'CA3') | ismember(cell_metrics.brainRegion,'CA3slm') | ismember(cell_metrics.brainRegion,'CA3slu') ...
            | ismember(cell_metrics.brainRegion,'CA3so') | ismember(cell_metrics.brainRegion,'CA3sp') | ismember(cell_metrics.brainRegion,'CA3sr') ...
            | ismember(cell_metrics.brainRegion,'CA2') | ismember(cell_metrics.brainRegion,'CA2slm') | ismember(cell_metrics.brainRegion,'CA2so') ...
            | ismember(cell_metrics.brainRegion,'CA2sp') | ismember(cell_metrics.brainRegion,'CA2sr') | ismember(cell_metrics.brainRegion,'DG') ...
            | ismember(cell_metrics.brainRegion,'HIP') | ismember(cell_metrics.brainRegion,'HPF');
        
        clear cellsID ccMedian ccZMedian ccMean ccZMean
        for jj = 1 : length(spikes.times)
            cellsID = indCell(indCell~=jj & cellsInRegion);

            ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
            ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG

            ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
            ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG
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

        brainRegionCCG.('HIP_medianCCG') = ccMedian;
        brainRegionCCG.('HIP_ZmedianCCG') = ccZMedian;
        brainRegionCCG.('HIP_meanCCG') = ccMean;
        brainRegionCCG.('HIP_ZmeanCCG') = ccZMean;
        brainRegionCCG.('HIP_ccgIndex') = ccgIndex;
        brainRegionsCcgIndex.('HIP') = ccgIndex;
        
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
if useDistinctShanks && length(session.extracellular.electrodeGroups.channels)>1 % if more than 1 shanks
   
    for ii =  1: length(session.extracellular.electrodeGroups.channels)
        cellsInShank = ismember(spikes.shankID,ii);
        
        clear cellsID ccMedian ccZMedian ccMean ccZMean
        for jj = 1 : length(spikes.times)
            cellsID = indCell(indCell~=jj & cellsInShank);
            
            if length(cellsID) == 1
                ccMedian(jj,:) = squeeze(allCcg(:,jj,cellsID));
                ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,cellsID))',[],2);
                
                ccMean(jj,:) = squeeze(allCcg(:,jj,cellsID));
                ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,cellsID))',[],2);
            else
                ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
                ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG

                ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
                ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2)); % zCCG
            end
            ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,cellsID)),2); %
            ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG

            ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,cellsID)),2); % zCCG
            ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,cellsID))',[],2),1); % zCCG
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

        shanksCCG.(['shank_' num2str(ii) '_medianCCG']) = ccMedian;
        shanksCCG.(['shank_' num2str(ii)  '_ZmedianCCG']) = ccZMedian;
        shanksCCG.(['shank_' num2str(ii)  '_meanCCG']) = ccMean;
        shanksCCG.(['shank_' num2str(ii)  '_ZmeanCCG']) = ccZMean;
        shanksCCG.(['shank_' num2str(ii)  '_ccgIndex']) = ccgIndex;
        shanksCCG.binSize = binSize;
        shanksCCG.winSize = winSize;
        shanksCCG.timestamps = t_ccg;
        shanksCCG.excludeIntervals = excludeIntervals;
        shanksCCG.winIndex = winIndex;
        shanksCCGIndex.(['shank_' num2str(ii)]) = ccgIndex;
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


% brainRegionCellTypeCCG = [];
% if useBrainRegions & useCellType & exist([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat'])
%     disp('Computing CCG by brain region and cell type');
%     session = loadSession;
%     load([basenameFromBasepath(basepath) '.cell_metrics.cellinfo.mat']);
%     
%     if isfield(session,'brainRegions')
%         efields = fieldnames(session.brainRegions);
%         is_pyr = strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell');
%         is_nw = strcmpi(cell_metrics.putativeCellType,'Narrow Interneuron');
%         is_ww = strcmpi(cell_metrics.putativeCellType,'Wide Interneuron');
%         is_int = strcmpi(cell_metrics.putativeCellType,'Narrow Interneuron') | strcmpi(cell_metrics.putativeCellType,'Wide Interneuron');
%         
%         for ii = 1:length(efields)
%             cellsInRegion = ismember(cell_metrics.brainRegion,efields{ii});
%             pyrCellsInRegion = cellsInRegion & is_pyr;
%             nwCellsInRegion = cellsInRegion & is_nw;
%             wwCellsInRegion = cellsInRegion & is_ww;
%             intCellsInRegion = cellsInRegion & is_int;
%             
%             clear cellsID ccMedian ccZMedian ccMean ccZMean
%             
%             % 1) pyrCellsID
%             for jj = 1:length(spikes.times)
%                pyrCellsID = indCell(indCell ~= jj & pyrCellsInRegion); 
%                if length(pyrCellsID) == 1
%                     ccMedian(jj,:) = squeeze(allCcg(:,jj,pyrCellsID)); %
%                     ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2); % zCCG
% 
%                     ccMean(jj,:) = squeeze(allCcg(:,jj,pyrCellsID)); % zCCG
%                     ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2); % zCCG
%                 else
%                     ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,pyrCellsID)),2); %
%                     ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2)); % zCCG
% 
%                     ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,pyrCellsID)),2); % zCCG
%                     ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2)); % zCCG
%                 end
%                 
%                 ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,pyrCellsID)),2); %
%                 ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2),1); % zCCG
% 
%                 ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,pyrCellsID)),2); % zCCG
%                 ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,pyrCellsID))',[],2),1); % zCCG
%                 
%                 if interp0
%                     artifactSamples = find(t_ccg == 0);
%                     x_axis = 1:length(t_ccg);
%                     x_axis(artifactSamples) = [];
%                     for jj = 1:size(ccMedian,1)
%                         ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
%                         ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
%                         ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
%                         ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
%                     end
%                 end
%                 
%                 win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
%                 ccgIndex = median(ccZMedian(:,win),2);
% 
%                 brainRegionCellTypeCCG.([efields{ii} '_pyr_medianCCG']) = ccMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_pyr_ZmedianCCG']) = ccZMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_pyr_meanCCG']) = ccMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_pyr_ZmeanCCG']) = ccZMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_pyr_ccgIndex']) = ccgIndex;
%                 brainRegionCellTypeCCG.binSize = binSize;
%                 brainRegionCellTypeCCG.winSize = winSize;
%                 brainRegionCellTypeCCG.timestamps = t_ccg;
%                 brainRegionCellTypeCCG.excludeIntervals = excludeIntervals;
%                 brainRegionCellTypeCCG.winIndex = winIndex;
%                 brainRegionCellTypeCcgIndex_pyr.(efields{ii}) = ccgIndex;
%                 
%                 
%             end
%             clear ccMedian ccZMedian ccMean ccZMean
%             
%             
%            % 2) nwCellsID            
%            for jj = 1:length(spikes.times)
%                
%                nwCellsID = indCell(indCell ~= jj & nwCellsInRegion);
%                
%                if length(nwCellsID) == 1
%                     ccMedian(jj,:) = squeeze(allCcg(:,jj,nwCellsID)); %
%                     ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2); % zCCG
% 
%                     ccMean(jj,:) = squeeze(allCcg(:,jj,nwCellsID)); % zCCG
%                     ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2); % zCCG
%                 else
%                     ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,nwCellsID)),2); %
%                     ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2)); % zCCG
% 
%                     ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,nwCellsID)),2); % zCCG
%                     ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2)); % zCCG
%                 end
%                 
%                 ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,nwCellsID)),2); %
%                 ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2),1); % zCCG
% 
%                 ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,nwCellsID)),2); % zCCG
%                 ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,nwCellsID))',[],2),1); % zCCG
%                 
%                 if interp0
%                     artifactSamples = find(t_ccg == 0);
%                     x_axis = 1:length(t_ccg);
%                     x_axis(artifactSamples) = [];
%                     for jj = 1:size(ccMedian,1)
%                         ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
%                         ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
%                         ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
%                         ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
%                     end
%                 end
%                 
%                 win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
%                 ccgIndex = median(ccZMedian(:,win),2);
% 
%                 brainRegionCellTypeCCG.([efields{ii} '_nw_medianCCG']) = ccMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_nw_ZmedianCCG']) = ccZMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_nw_meanCCG']) = ccMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_nw_ZmeanCCG']) = ccZMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_nw_ccgIndex']) = ccgIndex;
%                 brainRegionCellTypeCCG.binSize = binSize;
%                 brainRegionCellTypeCCG.winSize = winSize;
%                 brainRegionCellTypeCCG.timestamps = t_ccg;
%                 brainRegionCellTypeCCG.excludeIntervals = excludeIntervals;
%                 brainRegionCellTypeCCG.winIndex = winIndex;
%                 brainRegionCellTypeCcgIndex_nw.(efields{ii}) = ccgIndex;
%                 
%            end
%            clear ccMedian ccZMedian ccMean ccZMean
%                
%            % 3) wwCellsID    
%            for jj = 1:length(spikes.times)
%                
%                wwCellsID = indCell(indCell ~= jj & wwCellsInRegion);
%                
%                if length(wwCellsID) == 1
%                     ccMedian(jj,:) = squeeze(allCcg(:,jj,wwCellsID)); %
%                     ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2); % zCCG
% 
%                     ccMean(jj,:) = squeeze(allCcg(:,jj,wwCellsID)); % zCCG
%                     ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2); % zCCG
%                 else
%                     ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,wwCellsID)),2); %
%                     ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2)); % zCCG
% 
%                     ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,wwCellsID)),2); % zCCG
%                     ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2)); % zCCG
%                 end
%                 
%                 ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,wwCellsID)),2); %
%                 ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2),1); % zCCG
% 
%                 ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,wwCellsID)),2); % zCCG
%                 ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,wwCellsID))',[],2),1); % zCCG
%                 
%                 if interp0
%                     artifactSamples = find(t_ccg == 0);
%                     x_axis = 1:length(t_ccg);
%                     x_axis(artifactSamples) = [];
%                     for jj = 1:size(ccMedian,1)
%                         ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
%                         ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
%                         ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
%                         ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
%                     end
%                 end
%                 
%                 win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
%                 ccgIndex = median(ccZMedian(:,win),2);
% 
%                 brainRegionCellTypeCCG.([efields{ii} '_ww_medianCCG']) = ccMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_ww_ZmedianCCG']) = ccZMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_ww_meanCCG']) = ccMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_ww_ZmeanCCG']) = ccZMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_ww_ccgIndex']) = ccgIndex;
%                 brainRegionCellTypeCCG.binSize = binSize;
%                 brainRegionCellTypeCCG.winSize = winSize;
%                 brainRegionCellTypeCCG.timestamps = t_ccg;
%                 brainRegionCellTypeCCG.excludeIntervals = excludeIntervals;
%                 brainRegionCellTypeCCG.winIndex = winIndex;
%                 brainRegionCellTypeCcgIndex_ww.(efields{ii}) = ccgIndex;
%                 
%            end
%            clear ccMedian ccZMedian ccMean ccZMean
%                
%            % 4) intCellsID
%            for jj = 1:length(spikes.times)    
%                
%                intCellsID = indCell(indCell ~= jj & intCellsInRegion);
%                if length(intCellsID) == 1
%                     ccMedian(jj,:) = squeeze(allCcg(:,jj,intCellsID)); %
%                     ccZMedian(jj,:) = zscore(squeeze(allCcg(:,jj,intCellsID))',[],2); % zCCG
% 
%                     ccMean(jj,:) = squeeze(allCcg(:,jj,intCellsID)); % zCCG
%                     ccZMean(jj,:) = zscore(squeeze(allCcg(:,jj,intCellsID))',[],2); % zCCG
%                 else
%                     ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,intCellsID)),2); %
%                     ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,intCellsID))',[],2)); % zCCG
% 
%                     ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,intCellsID)),2); % zCCG
%                     ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,intCellsID))',[],2)); % zCCG
%                 end
%                 
%                 ccMedian(jj,:) = nanmedian(squeeze(allCcg(:,jj,intCellsID)),2); %
%                 ccZMedian(jj,:) = nanmedian(zscore(squeeze(allCcg(:,jj,intCellsID))',[],2),1); % zCCG
% 
%                 ccMean(jj,:) = nanmean(squeeze(allCcg(:,jj,intCellsID)),2); % zCCG
%                 ccZMean(jj,:) = nanmean(zscore(squeeze(allCcg(:,jj,intCellsID))',[],2),1); % zCCG
%                 
%                 if interp0
%                     artifactSamples = find(t_ccg == 0);
%                     x_axis = 1:length(t_ccg);
%                     x_axis(artifactSamples) = [];
%                     for jj = 1:size(ccMedian,1)
%                         ccMedian(jj,artifactSamples) = interp1(x_axis,ccMedian(jj,x_axis),artifactSamples);
%                         ccZMedian(jj,artifactSamples) = interp1(x_axis,ccZMedian(jj,x_axis),artifactSamples);
%                         ccMean(jj,artifactSamples) = interp1(x_axis,ccMean(jj,x_axis),artifactSamples);
%                         ccZMean(jj,artifactSamples) = interp1(x_axis,ccZMean(jj,x_axis),artifactSamples);
%                     end
%                 end
%                 
%                 win = t_ccg >= winIndex(1) & t_ccg <= winIndex(2);
%                 ccgIndex = median(ccZMedian(:,win),2);
% 
%                 brainRegionCellTypeCCG.([efields{ii} '_int_medianCCG']) = ccMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_int_ZmedianCCG']) = ccZMedian;
%                 brainRegionCellTypeCCG.([efields{ii} '_int_meanCCG']) = ccMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_int_ZmeanCCG']) = ccZMean;
%                 brainRegionCellTypeCCG.([efields{ii} '_int_ccgIndex']) = ccgIndex;
%                 brainRegionCellTypeCCG.binSize = binSize;
%                 brainRegionCellTypeCCG.winSize = winSize;
%                 brainRegionCellTypeCCG.timestamps = t_ccg;
%                 brainRegionCellTypeCCG.excludeIntervals = excludeIntervals;
%                 brainRegionCellTypeCCG.winIndex = winIndex;
%                 brainRegionCellTypeCcgIndex_int.(efields{ii}) = ccgIndex;
%            end
%            clear ccMedian ccZMedian ccMean ccZMean
%         end
%         averageCCG.brainRegionCellTypeCCG = brainRegionCellTypeCCG;
%     end
% end

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
        saveas(gcf,['SummaryFigures\grandCCGAverage.png']);
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
            saveas(gcf,['SummaryFigures\CCGAvgPerRegion.png']); 
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
            saveas(gcf,['SummaryFigures\CCGAvgPerShank.png']);
        end      
    end
    
end
close all;
cd(prevPath);
end