function [] = plotSummary_Drug(varargin)
%        [] = getSummaryPerCells(varargin)
%
% Display summary plots for a given cell in a single figures. By default
% goes over all cells of a given session
%
% INPUTS
% <Optional>
% 'basepath'            - Default, pwd
% 'UID'                 - Unique identifier for each neuron in a recording (see
%                           loadSpikes). If not provided, runs all cells.
% 'saveFigure'          - Default, true (in '/SummaryFigures/SummaryPerCell') 
% 'checkUnits'          - If true, ask which cells ID should be
%                           discarted...
% 'use_deltaThetaEpochs'- Default, true.
%
%% Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'UID',[], @isnumeric);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'checkUnits',true, @islogical);
addParameter(p,'use_deltaThetaEpochs',true, @islogical);
addParameter(p,'excludePlot',[]);
addParameter(p,'saveAs',[],@ischar);

parse(p,varargin{:})

basepath = p.Results.basepath;
UID = p.Results.UID;
saveFigure = p.Results.saveFigure;
checkUnits = p.Results.checkUnits;
use_deltaThetaEpochs = p.Results.use_deltaThetaEpochs;
excludePlot = p.Results.excludePlot;
saveAs = p.Results.saveAs;

for ii = 1:length(excludePlot)
    excludePlot{ii} = num2str(excludePlot{ii});
end
if length(excludePlot) == 0
    excludePlot = num2str(excludePlot);
end
excludePlot = lower(excludePlot);


% dealing with inputs 
prevPath = pwd;
cd(basepath);

if isempty(UID)
    spikes = loadSpikes;
    UID = spikes.UID;
end

% collecting pieces
targetFile = dir('*.cell_metrics.cellinfo.mat'); load(targetFile.name);
all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');

% Cell metrics Baseline
try
    targetFile = dir('*.cell_metrics_Drug.cellinfo.mat'); load(targetFile.name);
catch
    
end

% waveforms
all_waveforms = zscore(reshape([cell_metrics.waveforms.filt{:}],...
    [length(cell_metrics.waveforms.time{1}) length(cell_metrics.waveforms.filt)]));

pyr_color = [1 .7 .7];
nw_color = [.7 .7 1];
ww_color = [.7 1 1];

pyr_color_dark = pyr_color/4;
nw_color_dark = nw_color/4;
ww_color_dark = ww_color/4;

% acg
acg_time = [-50 : 0.5 : 50];
acg = cell_metrics.acg.narrow;
acg = acg./sum(acg);

% acg PeakTime
targetFile = dir('*.ACGPeak_Drug.cellinfo.mat'); 
acgPeak = importdata(targetFile.name);

targetFile = dir('*.thetaEpochs.states.mat'); load(targetFile.name);

% firing rate
spikemat = bz_SpktToSpkmat(loadSpikes,'dt',10,'units','rate');
if use_deltaThetaEpochs
    states_rate =   [cell_metrics.firingRate' cell_metrics.firingRate_QWake_noRipples_ThDt'  cell_metrics.firingRate_WAKEtheta_ThDt' cell_metrics.firingRate_NREM_noRipples_ThDt' cell_metrics.firingRate_REMtheta_ThDt'];
else
    states_rate =   [cell_metrics.firingRate' cell_metrics.firingRate_WAKEnontheta'          cell_metrics.firingRate_WAKEtheta'      cell_metrics.firingRate_NREMstate'           cell_metrics.firingRate_REMstate'];
end
run_quiet_index = (states_rate(:,3) - states_rate(:,2))./(states_rate(:,3) + states_rate(:,2));
rem_nrem_index = (states_rate(:,5) - states_rate(:,4))./(states_rate(:,5) + states_rate(:,4));

run_quiet_index = run_quiet_index/std(run_quiet_index);
rem_nrem_index = rem_nrem_index/std(rem_nrem_index);
cv = cell_metrics.firingRateCV; cv = cv/std(cv);

% avg CCG
targetFile = dir('*.averageCCG_Drug.cellinfo.mat'); load(targetFile.name);

% Ripples Responses Baseline
targetFile = dir('*.ripples_Drug_psth.cellinfo.mat'); ripplesResponses = importdata(targetFile.name);
targetFile = dir('*.ripple_120-200_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.ripples_Drug.events.mat'); load(targetFile.name);

% theta and gamma/s
targetFile = dir('*.theta_6-12_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.lgamma_20-60_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.hgamma_60-100_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

% Theta for REM and RUN
targetFile = dir('*.thetaRun_6-12_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.thetaREM_6-12_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

% speed (normalizing by firing rate)
targetFile = dir('*.speedCorrs_Drug.cellinfo.mat'); load(targetFile.name);
speedCorr = speedCorrs;
if ~isempty(speedCorr)
    for ii = 1:size(speedCorr.speedVals,1)
        speedVals(ii,:) = mean(speedCorr.speedVals(ii,:,:),3)/cell_metrics.firingRate(ii);
    end
end
    
% spatial modulation
targetFile = dir('*.spatialModulation.cellinfo.mat'); 
if ~isempty(targetFile)
    load(targetFile.name);
else 
    spatialModulation = [];
end

% behavioural events
targetFile = dir('*.behavior.cellinfo.mat');
if ~isempty(targetFile)
    load(targetFile.name);
else
    behavior = [];
end

session = loadSession;

if isempty(UID)
    showTagCells = false;
    UID = 1;
end
    

figure('units','normalized','outerposition',[0 0 1 1])
% waveform
subplot(6,6,1)
hold on
if length(find(all_nw) == 1) == 1
    plot(cell_metrics.waveforms.time{1}, all_waveforms(:,all_nw),'color',nw_color);
else
    plotFill(cell_metrics.waveforms.time{1}, all_waveforms(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.9);
end
if length(find(all_pyr) == 1) == 1
    plot(cell_metrics.waveforms.time{1}, all_waveforms(:,all_pyr),'color',pyr_color)
else
    plotFill(cell_metrics.waveforms.time{1}, all_waveforms(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.9);
end
if length(find(all_ww) == 1) == 1
    plot(cell_metrics.waveforms.time{1},all_waveforms(:,all_ww),'color',ww_color);
else
    plotFill(cell_metrics.waveforms.time{1}, all_waveforms(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.9);
end

scatter(cell_metrics.troughToPeak(all_pyr),rand(length(find(all_pyr)),1)/10 + 2,20,pyr_color,'filled');
scatter(cell_metrics.troughToPeak(all_nw),rand(length(find(all_nw)),1)/10 + 2.2,20,nw_color,'filled');
scatter(cell_metrics.troughToPeak(all_ww),rand(length(find(all_ww)),1)/10 + 2.4,20,ww_color,'filled');

axis tight; xlabel('ms'); ylabel('Waveform amp (SD)');

% acg
subplot(6,6,2)
hold on
if length(find(all_pyr) == 1) == 1
    plot(acg_time, acg(:,all_pyr),'color',pyr_color);
else
    plotFill(acg_time, acg(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.7);
end
if length(find(all_nw == 1)) == 1
    plot(acg_time, acg(:,all_nw),'color',nw_color);
else
    plotFill(acg_time, acg(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.7);
end
if length(find(all_ww) == 1) == 1
    plot(acg_time,acg(:,all_ww),'color',ww_color);
else
    plotFill(acg_time, acg(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.7);
end
scatter(rand(length(find(all_pyr)),1)*2 + 52, cell_metrics.acg_tau_rise(all_pyr)/1000,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)*2 + 55, cell_metrics.acg_tau_rise(all_nw)/1000,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)*2 + 58, cell_metrics.acg_tau_rise(all_ww)/1000,20,ww_color,'filled');

axis tight; xlabel('ms'); ylabel('ACG (prob)');

% ACG Peak
subplot(6,6,3)
hold on
if length(find(all_pyr) == 1) == 1
    plot(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_pyr),'color',pyr_color);
else
    plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.7);
end
if length(find(all_nw == 1)) == 1
    plot(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_nw),'color',nw_color);
else
    plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.7);
end
if length(find(all_ww) == 1) == 1
    plot(acgPeak.acg_time,acgPeak.acg_smoothed_norm(:,all_ww),'color',ww_color);
else
    plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.7);
end

scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_pyr)),rand(length(find(all_pyr)),1)/100 + 0.03,20,pyr_color,'filled');
scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_nw)),rand(length(find(all_nw)),1)/100 + 0.032,20,nw_color,'filled');
scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_ww)),rand(length(find(all_ww)),1)/100 + 0.036,20,ww_color,'filled');

XTick = [-2 -1 0 1];
set(gca,'XTick',XTick);
XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
set(gca,'XTickLabel',XTickLabels);
ylabel('logACG (prob)'); xlabel('Time(s)');
axis tight;

% cell position
subplot(6,6,4)
hold on
scatter(cell_metrics.general.chanCoords.x+100, cell_metrics.general.chanCoords.y,10,[.9 .9 .9],"filled");
scatter(cell_metrics.trilat_x(all_pyr)+100, cell_metrics.trilat_y(all_pyr),30,pyr_color,'x','LineWidth',2);
scatter(cell_metrics.trilat_x(all_nw)+100, cell_metrics.trilat_y(all_nw),30,nw_color,'x','LineWidth',2);
scatter(cell_metrics.trilat_x(all_ww)+100, cell_metrics.trilat_y(all_ww),30,ww_color,'x','LineWidth',2);
ax = axis;
ylabel('depth (\mum)'); xlabel('\mum')
edges = [min(cell_metrics.general.chanCoords.y)-50:40:max(cell_metrics.general.chanCoords.y)+50];
xHist = smooth(histcounts(cell_metrics.trilat_y(all_pyr), edges),1);
xHist = xHist/max(xHist); xHist = xHist *  300;
centers = edges(1:end-1) + mean(diff(edges))/2;
patch([0 xHist' 0],centers([1 1:end 1]),pyr_color,'EdgeColor','none','FaceAlpha',.5);
xHist = smooth(histcounts(cell_metrics.trilat_y(all_nw), edges),1);
xHist = xHist/max(xHist); xHist = xHist *  300;
patch([0 xHist' 0],centers([1 1:end 1]),nw_color,'EdgeColor','none','FaceAlpha',.5);
xHist = smooth(histcounts(cell_metrics.trilat_y(all_ww), edges),1);
xHist = xHist/max(xHist); xHist = xHist *  300;
patch([0 xHist' 0],centers([1 1:end 1]),ww_color,'EdgeColor','none','FaceAlpha',.5);
xlim([min(cell_metrics.general.chanCoords.x)+50 max(cell_metrics.general.chanCoords.x)+150]);
ylim([min(cell_metrics.general.chanCoords.y)-100 max(cell_metrics.general.chanCoords.y)+100]);

% Brain Region
subplot(6,6,5)
pie(categorical(cell_metrics.brainRegion))

% firing rate
subplot(6,6,6)
hold on
if length(find(all_pyr) == 1) == 1
    plot(spikemat.timestamps/60, spikemat.data(:,all_pyr),'color',pyr_color)
else
    plotFill(spikemat.timestamps/60, spikemat.data(:,all_pyr),'color',pyr_color,'yscale','log','style','filled','smoothOpt',10);
end
if length(find(all_nw == 1)) == 1
    plot(spikemat.timestamps/60, spikemat.data(:,all_nw),'color',nw_color);
else
    plotFill(spikemat.timestamps/60, spikemat.data(:,all_nw),'color',nw_color,'yscale','log','style','filled','smoothOpt',10);
end
if length(find(all_ww) == 1) == 1
    plot(spikemat.timestamps/60,spikemat.data(:,all_ww),'color',ww_color);
else
    plotFill(spikemat.timestamps/60, spikemat.data(:,all_ww),'color',ww_color,'yscale','log','style','filled','smoothOpt',10);
end

title(['Stability:', num2str(round(mean(cell_metrics.firingRateInstability),2)) char(177)...
    num2str(round(std(cell_metrics.firingRateInstability),2))],'FontWeight','normal');

ylim([.1 50]); xlim([0 spikemat.timestamps(end)/60]);
xlabel('Time (min)'); ylabel('Rate (Hz)');

% firing rate states
subplot(6,6,7)
hold on
if length(find(all_pyr) == 1) == 1
    plot(1:5,states_rate(all_pyr,:),'color',pyr_color)
else
    plotFill(1:5,states_rate(all_pyr,:),'color',pyr_color,'yscale','log','style','filled');
end
if length(find(all_nw == 1)) == 1
    plot(1:5,states_rate(all_nw,:),'color',nw_color);
else
    plotFill(1:5,states_rate(all_nw,:),'color',nw_color,'yscale','log','style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(1:5,states_rate(all_ww,:),'color',ww_color);
else
    plotFill(1:5,states_rate(all_ww,:),'color',ww_color,'yscale','log','style','filled');
end

set(gca,'XTick',[1:5],'XTickLabel',{'All','Wake','Run','NREM','REM'},'XTickLabelRotation',45);
xlim([.8 5.2]); ylabel('Rate (Hz)');

% States index
subplot(6,6,8)
hold on
scatter(rand(length(find(all_pyr)),1)/5 + .9, run_quiet_index(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 1.1, run_quiet_index(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 1.3, run_quiet_index(all_ww),20,ww_color,'filled');

scatter(rand(length(find(all_pyr)),1)/5 + 1.9, rem_nrem_index(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 2.1, rem_nrem_index(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 2.3, rem_nrem_index(all_ww),20,ww_color,'filled');

scatter(rand(length(find(all_pyr)),1)/5 + 2.9, cv(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 3.1, cv(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 3.3, cv(all_ww),20,ww_color,'filled');

plot([.8 3.4],[0 0],'-k');
set(gca,'TickDir','out','XTick',[1:3],'XTickLabel',{'run-quiet','rem-nrem','cv'},'XTickLabelRotation',45);
ylabel('Index (SD)'); xlim([.8 3.4]);

% avCCG
subplot(6,6,9)
hold on
if length(find(all_nw == 1)) == 1
    plot(averageCCG.timestamps,averageCCG.ZmedianCCG(all_nw,:),'color',nw_color);
else
    plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_nw,:),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(averageCCG.timestamps,averageCCG.ZmedianCCG(all_ww,:),'color',ww_color);
else
    plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_ww,:),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(averageCCG.timestamps,averageCCG.ZmedianCCG(all_pyr,:),'color',pyr_color)
else
    plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_pyr,:),'color',pyr_color,'style','filled');
end

scatter(rand(length(find(all_pyr)),1)/10 + .35, averageCCG.ccgIndex(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/10 + .4, averageCCG.ccgIndex(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/10 + .45, averageCCG.ccgIndex(all_ww),20,ww_color,'filled');

axis tight
xlabel('Time (s)'); ylabel('Rate (SD)');

% avCCG per region
subplot(6,6,10)
hold on
b = bar(averageCCG.brainRegionCCG.listOfRegionsID, nanmean(averageCCG.brainRegionCCG.absCcgIndexRegion));
b.FaceColor = 'flat';
b.FaceAlpha = 0.5;
b.EdgeColor = 'none';
b.CData(:,:) = hsv(length(averageCCG.brainRegionCCG.listOfRegionsID));
ylabel('CCG Index');
set(gca,'TickDir','out','XTick',averageCCG.brainRegionCCG.listOfRegionsID,'XTickLabel',averageCCG.brainRegionCCG.listOfRegions,'XTickLabelRotation',45);                  

% ripples 
subplot(6,6,11)
hold on
t_win = ripplesResponses.timestamps > -0.25 & ripplesResponses.timestamps < 0.25;
if length(find(all_nw == 1)) == 1
    plot(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_nw,t_win),'color',nw_color);
else
    plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveSmooth(all_ww,t_win),'color',ww_color);
else
    plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color)
else
    plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled','faceAlpha',.9);
end

scatter(rand(length(find(all_pyr)),1)/20 + .25, ripplesResponses.rateZDuringPulse(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/20 + .28, ripplesResponses.rateZDuringPulse(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/20 + .31, ripplesResponses.rateZDuringPulse(all_ww),20,ww_color,'filled');

title(['Pyr: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_pyr)),1)) ...
', NW: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_nw)),1))...
', WW: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_ww)),1)),' SD'],'FontWeight','normal');

plot(ripples.rippleStats.maps.timestamps, zscore(mean(ripples.rippleStats.maps.ripples_raw)) +20,'-k')
axis tight
xlabel('Ripple center (s)'); ylabel('Rate (SD)');

% ripple phase
subplot(6,6,12)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
if length(find(all_nw == 1)) == 1
    plot(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color);
else
    plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color);
else
    plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color)
else
    plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_pyr,:)','color',pyr_color,'style','filled');
end

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, rippleMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, rippleMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, rippleMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(rippleMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(rippleMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(rippleMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');

if isfield(session.analysisTags,'rippleChannel') & ~isempty(session.analysisTags.rippleChannel)
    for ii = 1:length(session.extracellular.electrodeGroups.channels)
          if ismember(session.analysisTags.rippleChannel,session.extracellular.electrodeGroups.channels{ii})
                ripple_shank = ii;  
          end
    end
else
    for ii = 1:length(session.extracellular.electrodeGroups.channels)
        if ismember(ripples.detectorinfo.detectionchannel,session.extracellular.electrodeGroups.channels{ii})
            ripple_shank = ii;
        end
    end     
end

flds = fields(session.brainRegions);
if isfield(session.analysisTags,'rippleChannel') & ~isempty(session.analysisTags.rippleChannel)
    for ii = 1:length(flds)    
        if ismember(session.analysisTags.rippleChannel,session.brainRegions.(flds{ii}).channels)
            ripple_region = flds{ii};
        end
    end
else
   for ii = 1:length(flds)    
        if ismember(ripples.detectorinfo.detectionchannel,session.brainRegions.(flds{ii}).channels)
            ripple_region = flds{ii};
        end
   end
end
    
title(['Ripple shank: ', num2str(ripple_shank), ' Ripple region: ', num2str(ripple_region) ],'FontWeight','normal');
axis tight
xlabel('Ripple phase (rad)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% theta
subplot(6,6,13)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
if length(find(all_nw == 1)) == 1
    plot(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color);
else
    plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color);
else
    plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color)
else
    plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:)','color',pyr_color,'style','filled');
end

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, thetaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, thetaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, thetaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(thetaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(thetaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(thetaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');

if isfield(session.analysisTags,'thetaChannel') & ~isempty(session.analysisTags.thetaChannel)
    for ii = 1:length(session.extracellular.electrodeGroups.channels)
        if ismember(session.analysisTags.thetaChannel,session.extracellular.electrodeGroups.channels{ii})
            theta_shank = ii;
        end
    end
else
    for ii = 1:length(session.extracellular.electrodeGroups.channels)
        if ismember(thetaEpochs.channel,session.extracellular.electrodeGroups.channels{ii})
            theta_shank = ii;
        end
    end
end
flds = fields(session.brainRegions);
if isfield(session.analysisTags,'thetaChannel') & ~isempty(session.analysisTags.thetaChannel)
    for ii = 1:length(flds)    
        if ismember(session.analysisTags.thetaChannel,session.brainRegions.(flds{ii}).channels)
            theta_region = flds{ii};
        end
    end
else
    for ii = 1:length(flds)    
        if ismember(thetaEpochs.channel,session.brainRegions.(flds{ii}).channels)
            theta_region = flds{ii};
        end
    end
end

title(['Theta shank: ', num2str(theta_shank), ' Theta region: ', num2str(theta_region) ],'FontWeight','normal');
axis tight
xlabel('Theta phase (rad)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});


% lgamma
subplot(6,6,14)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
if length(find(all_nw == 1)) == 1
    plot(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color);
else
    plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color);
else
    plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color)
else
    plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:)','color',pyr_color,'style','filled');
end

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, lgammaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, lgammaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, lgammaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(lgammaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(lgammaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(lgammaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');

axis tight
xlabel('Low gamma phase (rad)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% hgamma
subplot(6,6,15)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
if length(find(all_nw == 1)) == 1
    plot(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color);
else
    plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
end
if length(find(all_ww == 1)) == 1
    plot(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color);
else
    plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
end
if length(find(all_pyr) == 1) == 1
    plot(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color)
else
    plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:)','color',pyr_color,'style','filled');
end

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, hgammaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, hgammaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, hgammaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(hgammaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(hgammaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(hgammaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');

axis tight
xlabel('High gamma phase (rad)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% speed
subplot(6,6,16)
if ~isempty(speedCorr)
    if size(speedVals,1) ~= length(all_nw)
        speedVals = speedVals';
    end
    if length(find(all_nw == 1)) == 1
        plot(log10(speedCorr.prc_vals),speedVals(all_nw,:),'color',nw_color);
    else
        plotFill(log10(speedCorr.prc_vals),speedVals(all_nw,:),'color',nw_color,'style','filled');
    end
    if length(find(all_ww == 1)) == 1
        plot(log10(speedCorr.prc_vals),speedVals(all_ww,:),'color',ww_color);
    else
        plotFill(log10(speedCorr.prc_vals),speedVals(all_ww,:),'color',ww_color,'style','filled');
    end
    if length(find(all_pyr == 1)) ~= size(speedVals(all_pyr),2)
        plotFill(log10(speedCorr.prc_vals),speedVals(all_pyr,:),'color',pyr_color,'style','filled');
    else
        plotFill(log10(speedCorr.prc_vals),speedVals(all_pyr,:)','color',pyr_color,'style','filled');
    end
    scatter(rand(length(find(all_pyr)),1)/5 + 1.3, 1+speedCorr.speedScore(all_pyr)*2,20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/5 +  1.4, 1+speedCorr.speedScore(all_nw)*2,20,nw_color,'filled');
    scatter(rand(length(find(all_ww)),1)/5 + 1.5, 1+speedCorr.speedScore(all_ww)*2,20,ww_color,'filled');

    LogScale('x',10);
    axis tight
    xlabel('cm/s'); ylabel('speed rate/avg rate');
else
    axis off
    title('No speed data','FontWeight','normal');
end

% Behavior
if ~isempty(behavior)
    if isfield(behavior,'psth_entry') && isstruct(behavior.psth_entry)
        flds = fieldnames(behavior.psth_entry);
        for ii = 1:length(flds)
            for jj = 1:length(behavior.psth_entry.(flds{ii}))
                if ~isempty(behavior.psth_entry.(flds{ii}){jj})
                    h = get(gcf);
                    nextSubplot = length(h.Children) + 1;
                    subplot(6,6,nextSubplot)
                    hold on;
                    t_win = behavior.psth_entry.(flds{ii}){jj}.timestamps > -2 & behavior.psth_entry.(flds{ii}){jj}.timestamps < 2;
                    if length(find(all_nw == 1)) == 1 && length(find(all_nw == 1)) > 0
                        plot(behavior.psth_entry.(flds{ii}){jj}.timestamps(t_win),behavior.psth_entry.(flds{ii}){jj}.responsecurveZSmooth(all_nw,t_win),'color',nw_color);
                    elseif length(find(all_nw == 1)) > 1 
                        plotFill(behavior.psth_entry.(flds{ii}){jj}.timestamps(t_win),behavior.psth_entry.(flds{ii}){jj}.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
                    end
                    if length(find(all_ww == 1)) == 1 && length(find(all_ww == 1)) > 0
                        plot(behavior.psth_entry.(flds{ii}){jj}.timestamps(t_win),behavior.psth_entry.(flds{ii}){jj}.responsecurveZSmooth(all_ww,t_win),'color',ww_color);
                    elseif length(find(all_ww == 1)) > 1
                        plotFill(behavior.psth_entry.(flds{ii}){jj}.timestamps(t_win),behavior.psth_entry.(flds{ii}){jj}.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
                    end
                    plotFill(behavior.psth_entry.(flds{ii}){jj}.timestamps(t_win),behavior.psth_entry.(flds{ii}){jj}.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
                    axis tight
                    xlabel('Entry time (s)'); ylabel('Rate (SD)');
                    title([behavior.description{jj},' - ',flds{ii}],'FontWeight','normal');
                end
            end
        end
    end
    
    if isfield(behavior,'psth_exit') && isstruct(behavior.psth_exit)
        flds = fieldnames(behavior.psth_exit);
        for ii = 1:length(flds)
            for jj = 1:length(behavior.psth_exit.(flds{ii}))
                if ~isempty(behavior.psth_exit.(flds{ii}){jj})
                    h = get(gcf);
                    nextSubplot = length(h.Children) + 1;
                    subplot(6,6,nextSubplot)
                    hold on;
                    t_win = behavior.psth_exit.(flds{ii}){jj}.timestamps > -2 & behavior.psth_exit.(flds{ii}){jj}.timestamps < 2;
                    if length(find(all_nw == 1)) == 1
                        plot(behavior.psth_exit.(flds{ii}){jj}.timestamps(t_win),behavior.psth_exit.(flds{ii}){jj}.responsecurveZSmooth(all_nw,t_win),'color',nw_color);
                    else
                        plotFill(behavior.psth_exit.(flds{ii}){jj}.timestamps(t_win),behavior.psth_exit.(flds{ii}){jj}.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
                    end
                    if length(find(all_ww == 1)) == 1
                        plot(behavior.psth_exit.(flds{ii}){jj}.timestamps(t_win),behavior.psth_exit.(flds{ii}){jj}.responsecurveZSmooth(all_ww,t_win),'color',ww_color);
                    else
                        plotFill(behavior.psth_exit.(flds{ii}){jj}.timestamps(t_win),behavior.psth_exit.(flds{ii}){jj}.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
                    end
                    plotFill(behavior.psth_exit.(flds{ii}){jj}.timestamps(t_win),behavior.psth_exit.(flds{ii}){jj}.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
                    axis tight
                    xlabel('Exit time (s)'); ylabel('Rate (SD)');
                    title([behavior.description{jj},' - ',flds{ii}],'FontWeight','normal');
                end
            end
        end
    end
    
    if isfield(behavior,'psth_lReward') && isstruct(behavior.psth_lReward)
        h = get(gcf);
        nextSubPlot = length(h.Children) + 1;
        subplot(6,6,nextSubPlot)
        hold on;
        t_win = behavior.psth_lReward.timestamps > -2 & behavior.psth_lReward.timestamps < 2;
        if length(find(all_nw == 1)) == 1
            plot(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(all_nw,t_win),'color',nw_color);
        else
            plotFill(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
        end
        if length(find(all_ww == 1)) == 1
            plot(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(all_ww,t_win),'color',ww_color);
        else
            plotFill(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
        end
        plotFill(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
        axis tight
        xlabel('lReward time (s)'); ylabel('Rate (SD)');
    end
    
    if isfield(behavior,'psth_rReward') && isstruct(behavior.psth_rReward)
        h = get(gcf);
        nextSubPlot = length(h.Children) + 1;
        subplot(6,6,nextSubPlot)
        hold on;
        t_win = behavior.psth_rReward.timestamps > -2 & behavior.psth_rReward.timestamps < 2;
        if length(find(all_nw == 1)) == 1
            plot(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(all_nw,t_win),'color',nw_color);
        else
            plotFill(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
        end
        if length(find(all_ww == 1)) == 1
            plot(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(all_ww,t_win),'color',ww_color);
        else
            plotFill(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
        end
        plotFill(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
        axis tight
        xlabel('rReward time (s)'); ylabel('Rate (SD)');
    end
        
else
    h = get(gcf);
    nextSubPlot = length(h.Children) + 1;
    subplot(6,6,nextSubPlot)
    axis off
    title('No psth behaviour data','FontWeight','normal');
    
    h = get(gcf);
    nextSubPlot = length(h.Children) + 1;
    subplot(6,6,nextSubPlot)
    axis off
    title('No psth behaviour data','FontWeight','normal');
end

% spatial modulation
if ~isempty(spatialModulation) && ~any(ismember(excludePlot,{lower('spatialModulation')}))
    try
        h = get(gcf);
        nextSubPlot = length(h.Children) + 1;
        subplot(6,6,[nextSubplot nextSubplot+1])
        hold on

        imagesc_ranked(spatialModulation.map_1_timestamps, [1:length(find(all_pyr))], spatialModulation.map_1_rateMapsZ(all_pyr,:),[-3 3],...
            spatialModulation.PF_position_map_1(all_pyr,:));

        imagesc_ranked(spatialModulation.map_1_timestamps, [length(find(all_pyr)) + 5 (length(find(all_pyr))+ length(find(all_nw)) + 5)], spatialModulation.map_1_rateMapsZ(all_nw,:),[-3 3],...
            spatialModulation.PF_position_map_1(all_nw,:));

        imagesc_ranked(spatialModulation.map_1_timestamps,...
            [(length(find(all_pyr)) + length(find(all_nw)) + 10) (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)], spatialModulation.map_1_rateMapsZ(all_ww,:),[-3 3],...
            spatialModulation.PF_position_map_1(all_ww,:));
        xlim(spatialModulation.map_1_timestamps([1 end]));
        ylim([0 (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)]);
        xlabel('cm'); ylabel('Map 1 (-3 to 3 SD)');
        set(gca,'YTick',[0 (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)],'YTickLabel',[0 length(all_pyr)]);

        subplot(6,6,[23 24])
        hold on

        imagesc_ranked(spatialModulation.map_2_timestamps, [1:length(find(all_pyr))], spatialModulation.map_2_rateMapsZ(all_pyr,:),[-3 3],...
            spatialModulation.PF_position_map_2(all_pyr,:));

        imagesc_ranked(spatialModulation.map_2_timestamps, [length(find(all_pyr)) + 5 (length(find(all_pyr))+ length(find(all_nw)) + 5)], spatialModulation.map_2_rateMapsZ(all_nw,:),[-3 3],...
            spatialModulation.PF_position_map_1(all_nw,:));

        imagesc_ranked(spatialModulation.map_2_timestamps,...
            [(length(find(all_pyr)) + length(find(all_nw)) + 10) (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)], spatialModulation.map_2_rateMapsZ(all_ww,:),[-3 3],...
            spatialModulation.PF_position_map_1(all_ww,:));
        xlim(spatialModulation.map_2_timestamps([1 end]));
        ylim([0 (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)]);
        xlabel('cm'); ylabel('Map 2 (-3 to 3 SD)');
        set(gca,'YTick',[0 (length(find(all_pyr)) + length(find(all_nw)) + length(find(all_ww)) + 10)],'YTickLabel',[0 length(all_pyr)]);
        colormap jet
    catch
        subplot(6,6,[21 22])
        axis off
        title('No behaviour data','FontWeight','normal');

        subplot(6,6,[23 24])
        axis off
        title('No behaviour data','FontWeight','normal');
    end
else
    h = get(gcf);
    nextSubPlot = length(h.Children) + 1;
    subplot(6,6,[nextSubPlot])
    axis off
    title('No behaviour data','FontWeight','normal');
    
    h = get(gcf);
    nextSubPlot = length(h.Children) + 1;
    subplot(6,6,[nextSubPlot])
    axis off
    title('No behaviour data','FontWeight','normal');
end

if saveFigure
    saveas(gcf,['BaselinevsDrug\Summary_Drug.png']);
end

cd(prevPath);
end