function [spikesReturnPlot] = getSpikesReturnPlot(varargin)
% Computes Psth for timestamps entered as inputs.
% USAGE
%   [psth] = spikesPsth(<options>)
%
% INPUTS
%
% <OPTIONALS>
%   basepath - default pwd
%   spikes - buzcode spikes structure
%   doPlot - default, true
%   plot_caxis - if empty, uses max. Specifify otherwise (ex. [0 0.01])
% excludeStimulationPeriods
%               Default, true
% excludeIntervals 
%               2xN intervals to exlude
%
% OUTPUTS
%   spikesReturnPlot - struct
%
% Manu Valero 2022
%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'plot_caxis',[],@isnumeric);
addParameter(p,'plot_logxyaxis',[-3 2.5],@isnumeric);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'computeCellTypeAverages',true,@islogical);
addParameter(p,'saveMat',true,@islogical);

parse(p, varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
doPlot = p.Results.doPlot;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;
plot_caxis = p.Results.plot_caxis;
plot_logxyaxis = p.Results.plot_logxyaxis;
computeCellTypeAverages = p.Results.computeCellTypeAverages;
saveMat = p.Results.saveMat;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.spikesReturnPlot.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Return plot already computed! Loading file...');
    load(targetFile.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes();
end

if skipStimulationPeriods
    try
        optogenetic_responses = getOptogeneticResponse;
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

% compute return plot
for ii = 1:length(spikes.times)
    spikesReturnPlot.times{ii} = diff(spikes.times{ii}(1:end-1));
    spikesReturnPlot.times_plus1{ii} = diff(spikes.times{ii}(2:end));

    spikesReturnPlot.logTimes{ii} = log10(spikesReturnPlot.times{ii});
    spikesReturnPlot.logTimes_plus1{ii} = log10(spikesReturnPlot.times_plus1{ii});
    
    [N,c] = hist3([spikesReturnPlot.logTimes{ii} spikesReturnPlot.logTimes_plus1{ii}], 'Nbins',[100 100]);
    
    spikesReturnPlot.hist(ii,:,:) = N;
    spikesReturnPlot.hist_prob(ii,:,:) = N/sum(sum(N));
    spikesReturnPlot.centers = c{1};    
end

if computeCellTypeAverages
    try load([basenameFromBasepath(pwd) '.cell_metrics.cellinfo.mat']);
        cellTypes = unique(cell_metrics.putativeCellType);
        
        spikesReturnPlot.cellTypes = cell_metrics.putativeCellType;
        for ii = 1:length(cellTypes)
            spikesReturnPlot.(strrep(cellTypes{ii},' ','')).hist = ...
                squeeze(mean(spikesReturnPlot.hist(find(ismember(cell_metrics.putativeCellType,cellTypes{ii})),:,:),1));
            spikesReturnPlot.(strrep(cellTypes{ii},' ','')).hist_prob = ...
                squeeze(mean(spikesReturnPlot.hist_prob(find(ismember(cell_metrics.putativeCellType,cellTypes{ii})),:,:),1));
        end
    catch
        warning('Average returns maps per cell type not possible...');
    end
    
    cmap_pyr = makeColorMap([1 1 1],[1 .2 .2], [.7 0 0],50);
    cmap_nw = makeColorMap([1 1 1],[.2 .2 1], [0 0 .7],50);
    cmap_ww = makeColorMap([1 1 1],[.2 1 1], [0 .7 .7],50);
end

if doPlot
    % all cells
    figure;
    set(gcf,'Position',[200 -500 2500 1200]);
    for jj = 1:size(spikes.UID,2)
        % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        ax(jj) =subplot(7,ceil(size(spikes.UID,2)/7),jj);
        hold on
        contourf(spikesReturnPlot.centers,spikesReturnPlot.centers,squeeze(spikesReturnPlot.hist_prob(jj,:,:)),40,'LineColor','none');
        if ~isempty(plot_caxis)
            caxis(plot_caxis);
        else
            caxis([0 max(max(squeeze(spikesReturnPlot.hist_prob(jj,:,:))))]);
        end
        axis square
        ylim(plot_logxyaxis); xlim(plot_logxyaxis); 
        LogScale('xy',10);
        set(gca,'TickDir','out'); colormap(flip(gray));
        title(num2str(jj),'FontWeight','normal','FontSize',10);

        if jj == 1
            ylabel('ISI n+1 (s)');
            xlabel('ISI (s)');
        elseif jj == size(spikes.UID,2)
            xlabel('ISI (s)');
        else
            % set(gca,'YTick',[],'XTick',[]);
        end
    end
    saveas(gcf,['SummaryFigures\spikesReturnPlot.png']); 
    
    if isfield(spikesReturnPlot, 'cellTypes')
        for jj = 1:size(ax,2)
            switch spikesReturnPlot.cellTypes{jj}
                case 'Narrow Interneuron'
                    colormap(ax(jj),cmap_nw);
                case 'Wide Interneuron'
                    colormap(ax(jj),cmap_ww);
                case 'Pyramidal Cell'
                    colormap(ax(jj),cmap_pyr);
            end
        end
        saveas(gcf,['SummaryFigures\spikesReturnPlot.png']); 
            
        figure
        set(gcf,'Position',[100 100 900 400]);

        subplot(1,3,1);
        contourf(spikesReturnPlot.centers,spikesReturnPlot.centers,spikesReturnPlot.PyramidalCell.hist_prob,40,'LineColor','none');
        axis square
        ylim(plot_logxyaxis); xlim(plot_logxyaxis); 
        LogScale('xy',10);
        set(gca,'TickDir','out'); colormap(flip(gray));
        title('Pyramidal cells','FontWeight','normal','FontSize',10,'Color',[.8 .2 .2]);
        ylabel('ISI n+1 (s)');
        xlabel('ISI (s)');
        
        subplot(1,3,2);
        contourf(spikesReturnPlot.centers,spikesReturnPlot.centers,spikesReturnPlot.NarrowInterneuron.hist_prob,40,'LineColor','none');
        axis square
        ylim(plot_logxyaxis); xlim(plot_logxyaxis); 
        LogScale('xy',10);
        set(gca,'TickDir','out'); colormap(flip(gray));
        title('Narrow Waveform Int','FontWeight','normal','FontSize',10,'Color',[.2 .2 .8]);
        ylabel('ISI n+1 (s)');
        xlabel('ISI (s)');
        
        subplot(1,3,3);
        contourf(spikesReturnPlot.centers,spikesReturnPlot.centers,spikesReturnPlot.WideInterneuron.hist_prob,40,'LineColor','none');
        axis square
        ylim(plot_logxyaxis); xlim(plot_logxyaxis); 
        LogScale('xy',10);
        set(gca,'TickDir','out'); colormap(flip(gray));
        title('Wide Waveform Int','FontWeight','normal','FontSize',10,'Color',[.2 .8 .8]);
        ylabel('ISI n+1 (s)');
        xlabel('ISI (s)');
        saveas(gcf,['SummaryFigures\spikesReturnPlot_cellType.png']); 
            
    end
    
end

if saveMat
    disp('Saving results...');
    save([basenameFromBasepath(pwd) '.spikesReturnPlot.cellinfo.mat'],'spikesReturnPlot');
end

cd(prevPath);

end