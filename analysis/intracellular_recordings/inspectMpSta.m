
function [mpSta_restricted] = inspectMpSTA(varargin)
% inspectSTA(varargin)
%
% Inspection of membrante potential spike trigger-averages (STA)
%
% INPUTS
% <optional>
% basepath      Folder containing at least an .mpSta.intracellular.mat
%                   file (see getMpSta).
% mpSta         Matlab structure with mpSta results (see getMpSta).
% spikes        buzcode ripple structure (from bz_GetSpikes). If not provided, 
%                   it loads it from 'basepath' (if provided), or from current 
%                   folder (if not)
% restricSubthershold 
%               Show only subthreshold sweeps (true by default).
% membranePotentialInts
%               Restricto to specific membrane potential intervals. Default
%               [-Inf Inf]
% sr                (scalar) Sampling frequency in HZ (default loot at mpSta 
%               structure or 20000).
% winSTA            2 x 1 vector with STA time window (in ms), default [-10 100]
% saveSummary   Default true
% smooth        Default 10.
% guiSTA        Default false
%
% OUTPUT
% mpSta_restricted
%               mpSta structure as the one obtained by geMpSta but
%                   restricted according to the function inputs.
%  
%   Manu Valero 2018
%   Updated on 2020 for buzcode compatibility.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepah',pwd,@isnumeric);
addParameter(p,'mpSta',[],@isstruct);
addParameter(p,'spikes',[],@istruct);
addParameter(p,'sr',20000,@isscalar);
addParameter(p,'winSTA',[-100 100],@isnumeric);
addParameter(p,'restricSubthershold',true,@islogical);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'smoothFactor',10,@isscalar);
addParameter(p,'membranePotentialInts',[-Inf Inf], @isnumeric);
addParameter(p,'guiSTA',false,@islogical);

parse(p,varargin{:});
basepah = p.Results.basepah;
mpSta = p.Results.mpSta;
sr = p.Results.sr;
winSTA = p.Results.winSTA;
spikes = p.Results.spikes;
restricSubthershold = p.Results.restricSubthershold;
saveSummary = p.Results.saveSummary;
smoothFactor = p.Results.smoothFactor;
membranePotentialInts = p.Results.membranePotentialInts;
guiSTA = p.Results.guiSTA;

prevBasepath = pwd;
cd(basepah);

% dealing with inputs
if isempty(mpSta)
    disp('Loading STA data...');
    fileMpSta = dir('*mpSta.intracellular.mat');
    try load(fileMpSta.name,'mpSta');
    catch
        error('Intracellular file not found!');
    end
end
winSTA = winSTA/1000;

if isempty(spikes)
    spikes = loadSpikes;
end
xt = mpSta.timestamps;
idXt = find(xt >=winSTA(1) & xt <=winSTA(2));
intraCellNumbers = size(mpSta.data,2);
clusterCells = size(spikes.times,2);
[ccg, xccg]=CCG(spikes.times,[],'binSize',0.001,'duration',0.08,'norm','rate');

% RESTRICT MPSTA
mpSta_restricted = mpSta;
mpSta_restricted = rmfield(mpSta_restricted,'allTraces');
mpSta_restricted.restricSubthershold = restricSubthershold;
mpSta_restricted.membranePotentialInts = membranePotentialInts;
for ii = 1:size(intraCellNumbers)
    for jj = 1:size(spikes.times,2)
        if restricSubthershold 
            idx = find(mpSta.allTraces.subThre{ii}{jj});
        else
            idx = find(ones(size(mpSta.allTraces.subThre{ii}{jj})));
        end
        idMpInts = find(mpSta.allTraces.membranePotential{ii}{jj} > membranePotentialInts(1)...
            & mpSta.allTraces.membranePotential{ii}{jj} < membranePotentialInts(2));
        idx = intersect(idx, idMpInts);
        if numel(idx) > 10
            mpSta_restricted.data{ii}(:,jj) = mean(mpSta.allTraces.data{ii}{jj}(:,idx),2);
            mpSta_restricted.dataJittered{ii}(:,jj) = mean(mpSta.allTraces.data{ii}{jj}(:,idx),2)...
                        - mean(mpSta.allTraces.surr{ii}{jj}(idXt,:),2);
            mpSta_restricted.jitterStd{ii}(:,jj) = std(mpSta.allTraces.surr{ii}{jj}(idXt,:),[],2);
        else
            mpSta_restricted.data{ii}(:,jj) = nan(1,length(mpSta.timestamps));
            mpSta_restricted.dataJittered{ii}(:,jj) = nan(1,length(mpSta.timestamps));
            mpSta_restricted.jitterStd{ii}(:,jj) = nan(1,length(mpSta.timestamps));
        end
    end
end

% PLOT STA
if saveSummary
    mkdir('SummaryFigures'); % create folder
    % SAVE SUMMARY
    for ii = 1:size(intraCellNumbers)
        figure % STA
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.times,2)
            subplot(5,ceil(size(spikes.times,2)/5),jj); % autocorrelogram
            if size(mpSta.allTraces.data{ii}{jj},2) > 10
                if restricSubthershold 
                    idx = find(mpSta.allTraces.subThre{ii}{jj});
                else
                    idx = find(ones(size(mpSta.allTraces.subThre{ii}{jj})));
                end
                idMpInts = find(mpSta.allTraces.membranePotential{ii}{jj} > membranePotentialInts(1)...
                    & mpSta.allTraces.membranePotential{ii}{jj} < membranePotentialInts(2));
                idx = intersect(idx, idMpInts);
                
                if numel(idx) > 10
                    hold on;
                    plotFill(xt(idXt),mpSta.allTraces.surr{ii}{jj}(idXt,:),'color', [.8 .2 .2],'error','std','style','filled');
                    plot(xt(idXt), smooth(mean(mpSta.allTraces.data{ii}{jj}(:,idx),2), smoothFactor),'color',[.2 .2 .2],'LineWidth',1);
                    plot(xt(idXt), smooth(mean(mpSta.allTraces.data{ii}{jj}(:,idx),2)...
                        - mean(mpSta.allTraces.surr{ii}{jj}(idXt,:),2), smoothFactor),'color',[.2 .9 .2],'LineWidth',1);
                    xlim([min(xt(idXt)) max(xt(idXt))]); xlabel('ms'); ylabel('mV');
                    text(.05,.9,strcat(num2str(numel(idx)),' stas'),'Units','Normalized');
                    % title('STA','FontWeight','normal');
                    text(.05,.05,'STA','Units','Normalized','color',[.2 .2 .2]);
                    text(.3,.05,'Surr','Units','Normalized','color',[.8 .2 .2]);
                    text(.65,.05,'STA-Surr','Units','Normalized','color',[.2 .9 .2]);
                end
            end
        end
        saveas(gcf,['SummaryFigures\' 'mpStaCell_' num2str(membranePotentialInts(1)) 'to' num2str(membranePotentialInts(2)) '_'  num2str(ii) '.png']);
    
        figure % pre units
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.times,2)
            subplot(5,ceil(size(spikes.times,2)/5),jj); % autocorrelogram
            area(xccg,ccg(:,jj,jj),'LineStyle','none');
            xlim([min(xccg) max(xccg)]); xlabel('ms'); ylabel('Hz');
            text(.05,.9,strcat(num2str(numel(spikes.times{jj})),' spks, ',...
                num2str(round(numel(spikes.times{jj})/((spikes.times{jj}(end)-spikes.times{jj}(1))),2)),'Hz'),'Units','Normalized');
            title(strcat('Intra ','{ }', num2str(ii),'{ }', 'of','{ }', num2str(intraCellNumbers),...
                ', unit', '{ }', num2str(jj),'{ }', 'of','{ }', num2str(size(spikes.times,2))),'FontWeight','normal');
        end
        saveas(gcf,['SummaryFigures\' 'preUnits_mpStaCell_'  num2str(membranePotentialInts(1)) 'to' num2str(membranePotentialInts(2)) '_'  num2str(ii) '.png']);
        
        figure % all sweeps
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.times,2)
            subplot(5,ceil(size(spikes.times,2)/5),jj); % autocorrelogram
            if size(mpSta.allTraces.data{ii}{jj},2) > 10
                if restricSubthershold 
                    idx = find(mpSta.allTraces.subThre{ii}{jj});
                else
                    idx = find(ones(size(mpSta.allTraces.subThre{ii}{jj})));
                end
                idMpInts = find(mpSta.allTraces.membranePotential{ii}{jj} > membranePotentialInts(1)...
                    & mpSta.allTraces.membranePotential{ii}{jj} < membranePotentialInts(2));
                idx = intersect(idx, idMpInts);
                
                if numel(idx) > 10
                    [sortResting, idxResting] = sort(mpSta.allTraces.membranePotential{ii}{jj}(idx));
                    temp = mpSta.allTraces.data{ii}{jj}(idXt,idx);
                    imagesc(xt,1:length(idxResting),temp(:,idxResting)',[-2 2]);
                    xlim([min(xt(idXt)) max(xt(idXt))]); set(gca,'YDir','normal');
                    colormap jet
                    xlabel('ms'); ylabel('STAs (#)');
                end
            end
        end
        saveas(gcf,['SummaryFigures\' 'allSweeps_mpStaCell_' num2str(membranePotentialInts(1)) 'to' num2str(membranePotentialInts(2)) '_'  num2str(ii) '.png']);
    end
end

disp('Inspect intracell STA...');
if guiSTA
    figure
    ii = 1;
    while ii <= intraCellNumbers
        jj = 1;
        while jj <= size(spikes.times,2)
                disp('Plotting... ');

                clf;
                subplot(2,2,1)
                area(xccg,ccg(:,jj,jj),'LineStyle','none');
                xlim([min(xccg) max(xccg)]); xlabel('ms'); ylabel('Hz');
                text(.05,.9,strcat(num2str(numel(spikes.times{jj})),' spks, ',...
                    num2str(round(numel(spikes.times{jj})/((spikes.times{jj}(end)-spikes.times{jj}(1))),2)),'Hz'),'Units','Normalized');
                title(strcat('Intra ','{ }', num2str(ii),'{ }', 'of','{ }', num2str(intraCellNumbers),...
                    ', unit', '{ }', num2str(jj),'{ }', 'of','{ }', num2str(size(spikes.times,2))),'FontWeight','normal');

                if size(mpSta.allTraces.data{ii}{jj},2) > 10
                    if restricSubthershold 
                        idx = find(mpSta.allTraces.subThre{ii}{jj});
                    else
                        idx = find(ones(size(mpSta.allTraces.subThre{ii}{jj})));
                    end
                    if numel(idx) > 10
                        subplot(2,2,2)
                        hold on;
                        plotFill(xt(idXt),mpSta.allTraces.surr{ii}{jj}(idXt,:),'color', [.8 .2 .2],'error','std','style','filled');
                        plot(xt(idXt), smooth(mpSta.data{ii}(idXt,jj),smoothFactor),'color',[.2 .2 .2],'LineWidth',1);
                        plot(xt(idXt), smooth(mpSta.dataJittered{ii}(idXt,jj),smoothFactor),'color',[.2 .9 .2],'LineWidth',1);
                        xlim([min(xt(idXt)) max(xt(idXt))]); xlabel('ms'); ylabel('mV');
                        text(.05,.9,strcat(num2str(numel(idx)),' stas'),'Units','Normalized');
                        title('STA','FontWeight','normal');
                        text(.05,.05,'STA','Units','Normalized','color',[.2 .2 .2]);
                        text(.3,.05,'Surr','Units','Normalized','color',[.8 .2 .2]);
                        text(.65,.05,'STA-Surr','Units','Normalized','color',[.2 .9 .2]);

                        subplot(2,2,3)
                        [sortResting, idxResting] = sort(mpSta.allTraces.membranePotential{ii}{jj}(idx));
                        temp = mpSta.allTraces.data{ii}{jj}(idXt,idx);
                        imagesc(xt,1:length(idxResting),temp(:,idxResting)',[-2 2]);
                        xlim([min(xt(idXt)) max(xt(idXt))]); set(gca,'YDir','normal');
                        colormap jet
                        xlabel('ms'); ylabel('STAs (#)');
                        title('All STAs','FontWeight','normal');

                        subplot(2,2,4)
                        plot(sortResting,1:length(idxResting)); ylim([1 length(idxResting)]);
                        title('All STAs membrane potential','FontWeight','normal');
                        xlabel('mV'); ylabel('STAs (#)');
                    end
                end

                tog = input('Press ENTER or 6 to continue, 4 to go back: ');
                if isempty(tog) || tog == 6 
                    jj = jj + 1;
                elseif tog == 4
                    jj = jj - 1;
                else
                    warning('Input not recognized. ');
                end

        end
        fprintf('\nEnd of cell #%i.\n',ii);
        toh = input('Press ENTER or 6 to continue, 4 to go back: ');
        if isempty(tog) || tog == 6 
            ii = ii + 1;
        elseif tog == 4
            ii = ii - 1;
        else
            warning('Input not recognized. ');
        end
    end
end

cd(prevBasepath);
end