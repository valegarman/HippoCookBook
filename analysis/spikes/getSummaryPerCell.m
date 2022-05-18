function [] = getSummaryPerCell(varargin)
%        [] = getSummaryPerCells(varargin)
%
% Display summary plots for a given cell in a single figures. By default
% goes over all cells of a given session
%
% INPUTS
% <Optional>
% 'basepath'    - Default, pwd
% 'UID'         - Unique identifier for each neuron in a recording (see
%                   loadSpikes). If not provided, runs all cells.
% 'saveFigure'  - Default, true (in '/SummaryFigures/SummaryPerCell') 
% 'onlyOptoTag' - Runs code only in those cells with optogenetic responses.
%                   By default, true.
%
%% Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'UID',[], @isnumeric);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'onlyOptoTag',true, @islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
UID = p.Results.UID;
saveFigure = p.Results.saveFigure;
onlyOptoTag = p.Results.onlyOptoTag;

% dealing with inputs 
prevPath = pwd;
cd(basepath);

if onlyOptoTag
    optogenetic_responses = getOptogeneticResponse;
    UID = find(optogenetic_responses.threeWaysTest==1);
    clear optogenetic_responses
end

if isempty(UID)
    spikes = loadSpikes;
    UID = spikes.UID;
end

% collecting pieces
targetFile = dir('*.cell_metrics.cellinfo.mat'); load(targetFile.name);
all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
% waveforms
all_waveforms = zscore(reshape([cell_metrics.waveforms.filt{:}],...
    [length(cell_metrics.waveforms.time{1}) length(cell_metrics.waveforms.filt)]));
pyr_color = [1 .8 .8];
nw_color = [.8 .8 1];
cell_color = [0 0 0];
% acg
acg_time = [-50 : 0.5 : 50];
acg = cell_metrics.acg.narrow;
acg = acg./sum(acg);
% firing rate
spikemat = bz_SpktToSpkmat(loadSpikes,'dt',10,'units','rate');
states_rate = [cell_metrics.firingRate' cell_metrics.firingRate_WAKEnontheta' cell_metrics.firingRate_WAKEtheta' cell_metrics.firingRate_NREMstate' cell_metrics.firingRate_REMstate'];
run_quiet_index = (states_rate(:,3) - states_rate(:,2))./(states_rate(:,3) + states_rate(:,2));
rem_nrem_index = (states_rate(:,5) - states_rate(:,4))./(states_rate(:,5) + states_rate(:,4));

run_quiet_index = run_quiet_index/std(run_quiet_index);
rem_nrem_index = rem_nrem_index/std(rem_nrem_index);
cv = cell_metrics.firingRateCV; cv = cv/std(cv);

% avg CCG
[averageCCG] = getAverageCCG;

% theta and gamma/s
targetFile = dir('*.theta_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.lgamma_20-60.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.hgamma_60-100.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

for ii = 1:length(UID)
    figure;
    set(gcf,'Position',[50 1000 1400 800]);
    
    % waveform
    subplot(6,6,1)
    hold on
    plotFill(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,all_pyr),'style','filled','color',pyr_color);
    plotFill(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,all_nw),'style','filled','color',nw_color);
    plot(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,UID(ii)),'LineWidth',1.5,'color',cell_color);
    axis tight; xlabel('ms'); ylabel('Waveform amp (SD)');

    % acg
    subplot(6,6,2)
    hold on
    plotFill(acg_time, acg(:,all_nw),'style','filled','color',nw_color);
    plotFill(acg_time, acg(:,all_pyr),'style','filled','color',pyr_color);
    plot(acg_time, acg(:,UID(ii)),'LineWidth',1.5,'color',cell_color);
    axis tight; xlabel('ms'); ylabel('ACG (prob)');

    % cell position
    subplot(6,6,[3])
    hold on
    scatter(cell_metrics.general.chanCoords.x+100, cell_metrics.general.chanCoords.y,10,[.9 .9 .9],"filled");
    scatter(cell_metrics.trilat_x(all_pyr)+100, cell_metrics.trilat_y(all_pyr),30,pyr_color,'x','LineWidth',2);
    scatter(cell_metrics.trilat_x(all_nw)+100, cell_metrics.trilat_y(all_nw),30,nw_color,'x','LineWidth',2);
    scatter(cell_metrics.trilat_x(UID(ii))+100, cell_metrics.trilat_y(UID(ii)),50,cell_color,'x','LineWidth',2);
    ax = axis;
    ylabel('depth (\mum)'); xlabel('\mum')
    edges = [min(cell_metrics.general.chanCoords.y)-50:40:max(cell_metrics.general.chanCoords.y)+50];
    xHist = smooth(histcounts(cell_metrics.trilat_y(all_pyr), edges),1);
    xHist = xHist/max(xHist); xHist = xHist *  300;
    centers = edges(1:end-1) + mean(diff(edges))/2;
    patch([0 xHist' 0],centers([1 1:end 1]),pyr_color,'EdgeColor','none','FaceAlpha',.5);
    xHist = smooth(histcounts(cell_metrics.trilat_y(all_nw), edges),1);
    xHist = xHist/max(xHist); xHist = xHist *  300;
    centers = edges(1:end-1) + mean(diff(edges))/2;
    patch([0 xHist' 0],centers([1 1:end 1]),nw_color,'EdgeColor','none','FaceAlpha',.5);
    plot([0 300],[cell_metrics.trilat_y(UID(ii)) cell_metrics.trilat_y(UID(ii))],'color',cell_color,'LineWidth',1.5)
    title(cell_metrics.brainRegion(UID(ii)),'FontWeight','normal');
    
    % firing rate
    subplot(6,6,[4])
    hold on
    plotFill(spikemat.timestamps/60, spikemat.data(:,all_pyr),'color',pyr_color,'yscale','log','style','filled','smoothOpt',10);
    plotFill(spikemat.timestamps/60, spikemat.data(:,all_nw),'color',nw_color,'yscale','log','style','filled','smoothOpt',10);
    plot(spikemat.timestamps/60, smooth(spikemat.data(:,UID(ii)),10),'color',cell_color,'LineWidth',1);
    ylim([.1 50]); xlim([0 spikemat.timestamps(end)/60]);
    xlabel('Time (min)'); ylabel('Rate (Hz)');
    title(['Rate:' num2str(round(cell_metrics.firingRate(UID(ii)),2)),...
        'Hz, stability: ',  num2str(round(cell_metrics.firingRateInstability(UID(ii)),2)),...
        ' (Avg:', num2str(round(mean(cell_metrics.firingRateInstability),2)) char(177)...
        num2str(round(std(cell_metrics.firingRateInstability),2)),')'],'FontWeight','normal');

    % firing rate states
    subplot(6,6,[5])
    hold on
    plotFill(1:5,states_rate(all_pyr,:),'color',pyr_color,'yscale','log','style','filled');
    plotFill(1:5,states_rate(all_nw,:),'color',nw_color,'yscale','log','style','filled');
    plot(1:5,states_rate(UID(ii),:),'color',cell_color,'LineWidth',1);
    set(gca,'XTick',[1:5],'XTickLabel',{'All','Wake','Run','NREM','REM'},'XTickLabelRotation',45);

    subplot(6,6,[6])
    hold on
    scatter(rand(length(find(all_pyr)),1)/5 + .9, run_quiet_index(all_pyr),20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/5 + 1.1, run_quiet_index(all_nw),20,nw_color,'filled');
    scatter(1, run_quiet_index(UID(1)),40,cell_color,'filled');

    scatter(rand(length(find(all_pyr)),1)/5 + 1.9, rem_nrem_index(all_pyr),20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/5 + 2.1, rem_nrem_index(all_nw),20,nw_color,'filled');
    scatter(2, rem_nrem_index(UID(1)),40,cell_color,'filled');
    plot([.8 3.4],[0 0],'-');
    set(gca,'TickDir','out','XTick',[1:5],'XTickLabel',{'run-quiet','rem-nrem'},'XTickLabelRotation',45);
    ylabel('SD');

    scatter(rand(length(find(all_pyr)),1)/5 + 2.9, cv(all_pyr),20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/5 + 3.1, cv(all_nw),20,nw_color,'filled');
    scatter(3, cv(UID(1)),40,cell_color,'filled');
    set(gca,'TickDir','out','XTick',[1:3],'XTickLabel',{'run-quiet','rem-nrem','cv'},'XTickLabelRotation',45);
    ylabel('SD'); xlim([.8 3.4]);

    % avCCG
    subplot(6,6,7)
    hold on
    plotFill(averageCCG.timestamps,averageCCG.ZmeanCCG(all_nw,:),'color',nw_color,'style','filled');
    plotFill(averageCCG.timestamps,averageCCG.ZmeanCCG(all_pyr,:),'color',pyr_color,'style','filled');
    plot(averageCCG.timestamps,averageCCG.ZmeanCCG(UID(ii),:),'color',cell_color,'LineWidth',1.5);

    scatter(rand(length(find(all_pyr)),1)/10 + .35, averageCCG.ccgIndex(all_pyr),20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/10 + .4, averageCCG.ccgIndex(all_nw),20,nw_color,'filled');
    scatter(.375, run_quiet_index(UID(1)),40,cell_color,'filled');
    xlabel('Time (s)'); ylabel('SD');

    % theta... % complete with 
    subplot(6,6,8)
    hold on
    plotFill(thetaMod.timestamps,averageCCG.ZmeanCCG(all_nw,:),'color',nw_color,'style','filled');
    
    
    % gamma


    subplot(6,6,8)



    

    
    





    



end

cd(prevPath);
end
