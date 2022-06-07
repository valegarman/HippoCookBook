function [] = getSummaryPerSession(varargin)
%        [] = getSummaryPerSession(varargin)
%
% Display summary plots for all cells in a single figures. 
%
% INPUTS
% <Optional>
% 'basepath'            - Default, pwd
% 'UID'                 - Unique identifier for each neuron in a recording (see
%                           loadSpikes). If not provided, runs all cells.
% 'saveFigure'          - Default, true (in '/SummaryFigures/SummaryPerCell') 
% 'onlyOptoTag'         - Runs code only in those cells with optogenetic responses.
%                           By default, true.
% 'lightPulseDuration'  - In seconds, default 0.1.
%
%% Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'UID',[], @isnumeric);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'lightPulseDuration',0.1, @isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
UID = p.Results.UID;
saveFigure = p.Results.saveFigure;
lightPulseDuration = p.Results.lightPulseDuration;

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

% opto responses
optogenetic_responses = getOptogeneticResponse;
spikes = loadSpikes;
winSizePlot_opto = [-.1 .5];

% waveforms
all_waveforms = zscore(reshape([cell_metrics.waveforms.filt{:}],...
    [length(cell_metrics.waveforms.time{1}) length(cell_metrics.waveforms.filt)]));
pyr_color = [1 .7 .7];
nw_color = [.7 .7 1];
ww_color = [.7 1 1];
cell_color = [0 0 0];

pyr_color_dark = pyr_color/4;
nw_color_dark = nw_color/4;
ww_color_dark = ww_color/4;

% acg
acg_time = [-50 : 0.5 : 50];
acg = cell_metrics.acg.narrow;
acg = acg./sum(acg);

% acg PeakTime
targetFile = dir('*.ACGPeak.cellinfo.mat'); 
acgPeak = importdata(targetFile.name);

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

% ripples responses
targetFile = dir('*.ripples_psth.cellinfo.mat'); ripplesResponses = importdata(targetFile.name);
targetFile = dir('*.ripple_120-200.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
ripples = rippleMasterDetector;

% theta and gamma/s
targetFile = dir('*.theta_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.lgamma_20-60.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.hgamma_60-100.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

% speed
speedCorr = getSpeedCorr;
if ~isempty(speedCorr)
    for ii = 1:size(speedCorr.speedVals,2)
        speedVals(ii,:) = mean(speedCorr.speedVals(:,ii,:),3)/cell_metrics.firingRate(ii);
    end
end
    
% speedVals = bsxfun(@rdivide,mean(speedCorr.speedVals,3),cell_metrics.firingRate(:));

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

figure;
set(gcf,'Position',get(0,'screensize'));
% opto response
subplot(5,5,1)

% waveform
subplot(5,5,2)
hold on
plotFill(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.9);
plotFill(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.9);
plotFill(cell_metrics.waveforms.time{UID(ii)}, all_waveforms(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.9);

scatter(cell_metrics.troughToPeak(all_pyr),rand(length(find(all_pyr)),1)/10 + 2,20,pyr_color,'filled');
scatter(cell_metrics.troughToPeak(all_nw),rand(length(find(all_nw)),1)/10 + 2.2,20,nw_color,'filled');
scatter(cell_metrics.troughToPeak(all_ww),rand(length(find(all_ww)),1)/10 + 2.4,20,ww_color,'filled');
axis tight; xlabel('ms'); ylabel('Waveform amp (SD)');

% acg
subplot(5,5,3)
hold on
plotFill(acg_time, acg(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.7);
plotFill(acg_time, acg(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.7);
plotFill(acg_time, acg(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.7);
scatter(rand(length(find(all_pyr)),1)*2 + 52, cell_metrics.acg_tau_rise(all_pyr)/1000,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)*2 + 55, cell_metrics.acg_tau_rise(all_nw)/1000,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)*2 + 58, cell_metrics.acg_tau_rise(all_ww)/1000,20,ww_color,'filled');
axis tight; xlabel('ms'); ylabel('ACG (prob)');

% cell position
subplot(5,5,4)
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
title(cell_metrics.brainRegion(UID(ii)),'FontWeight','normal');
title([cell_metrics.brainRegion{UID(ii)} ', ch: '  num2str(cell_metrics.maxWaveformCh1(UID(ii)))],'FontWeight','normal');

% ACG Peak
subplot(5,5,5)
hold on
plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_pyr),'style','filled','color',pyr_color,'faceAlpha',0.7);
plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_nw),'style','filled','color',nw_color,'faceAlpha',0.7);
plotFill(acgPeak.acg_time, acgPeak.acg_smoothed_norm(:,all_ww),'style','filled','color',ww_color,'faceAlpha',0.7);

scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_pyr)),rand(length(find(all_pyr)),1)/100 + 0.04,20,pyr_color,'filled');
scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_nw)),rand(length(find(all_nw)),1)/100 + 0.042,20,nw_color,'filled');
scatter(acgPeak.acg_time(acgPeak.acgPeak_sample(all_ww)),rand(length(find(all_ww)),1)/100 + 0.046,20,ww_color,'filled');

% scatter(rand(length(find(all_pyr)),1)*0.1 + 1.02, acgPeak.acgPeak_sample(all_pyr)/1000,20,pyr_color,'filled');
% scatter(rand(length(find(all_nw)),1)*0.1 + 1.12, acgPeak.acgPeak_sample(all_nw)/1000,20,nw_color,'filled');
% scatter(rand(length(find(all_ww)),1)*0.1 + 1.22, acgPeak.acgPeak_sample(all_ww)/1000,20,ww_color,'filled');
% scatter(rand(length(find(UID(ii))),1)*0.1 + 1.32, acgPeak.acgPeak_sample(UID(ii))/1000,20,cell_color,'filled');

XTick = [-2 -1 0 1];
set(gca,'XTick',XTick);
XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
set(gca,'XTickLabel',XTickLabels);
ylabel('logACG (prob)'); xlabel('Time(s)');
axis tight;

% firing rate
subplot(5,5,6)
hold on
plotFill(spikemat.timestamps/60, spikemat.data(:,all_pyr),'color',pyr_color,'yscale','log','style','filled','smoothOpt',10);
plotFill(spikemat.timestamps/60, spikemat.data(:,all_nw),'color',nw_color,'yscale','log','style','filled','smoothOpt',10);
plotFill(spikemat.timestamps/60, spikemat.data(:,all_ww),'color',ww_color,'yscale','log','style','filled','smoothOpt',10);
ylim([.1 50]); xlim([0 spikemat.timestamps(end)/60]);
xlabel('Time (min)'); ylabel('Rate (Hz)');

% firing rate states
subplot(5,5,7)
hold on
plotFill(1:5,states_rate(all_pyr,:),'color',pyr_color,'yscale','log','style','filled');
plotFill(1:5,states_rate(all_nw,:),'color',nw_color,'yscale','log','style','filled');
plotFill(1:5,states_rate(all_ww,:),'color',ww_color,'yscale','log','style','filled');
set(gca,'XTick',[1:5],'XTickLabel',{'All','Wake','Run','NREM','REM'},'XTickLabelRotation',45);
xlim([.8 5.2]); ylabel('Rate (Hz)');

subplot(5,5,[8])
hold on
scatter(rand(length(find(all_pyr)),1)/5 + .9, run_quiet_index(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 1.1, run_quiet_index(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 1.3, run_quiet_index(all_ww),20,ww_color,'filled');

scatter(rand(length(find(all_pyr)),1)/5 + 1.9, rem_nrem_index(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 2.1, rem_nrem_index(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 2.3, rem_nrem_index(all_ww),20,ww_color,'filled');
plot([.8 3.4],[0 0],'-k');
set(gca,'TickDir','out','XTick',[1:5],'XTickLabel',{'run-quiet','rem-nrem'},'XTickLabelRotation',45);

scatter(rand(length(find(all_pyr)),1)/5 + 2.9, cv(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/5 + 3.1, cv(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/5 + 3.3, cv(all_ww),20,ww_color,'filled');
set(gca,'TickDir','out','XTick',[1:3],'XTickLabel',{'run-quiet','rem-nrem','cv'},'XTickLabelRotation',45);
ylabel('Index (SD)'); xlim([.8 3.4]);

% avCCG
subplot(5,5,9)
hold on
plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_nw,:),'color',nw_color,'style','filled');
plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_ww,:),'color',ww_color,'style','filled');
plotFill(averageCCG.timestamps,averageCCG.ZmedianCCG(all_pyr,:),'color',pyr_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/10 + .35, averageCCG.ccgIndex(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/10 + .4, averageCCG.ccgIndex(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/10 + .45, averageCCG.ccgIndex(all_ww),20,ww_color,'filled');
axis tight
xlabel('Time (s)'); ylabel('Rate (SD)');

% ripples
subplot(5,5,11)
hold on
t_win = ripplesResponses.timestamps > -0.25 & ripplesResponses.timestamps < 0.25;
plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled','faceAlpha',.9);
plotFill(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/20 + .25, ripplesResponses.rateZDuringPulse(all_pyr),20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/20 + .28, ripplesResponses.rateZDuringPulse(all_nw),20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/20 + .31, ripplesResponses.rateZDuringPulse(all_ww),20,ww_color,'filled');
plot(ripples.rippleStats.maps.timestamps, zscore(mean(ripples.rippleStats.maps.ripples_raw))+20,'-k')
axis tight
xlabel('Ripple center (s)'); ylabel('Rate (SD)');
title(['Pyr: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_pyr)),1)) ...
    ', NW: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_nw)),1))...
    ', WW: '  num2str(round(mean(ripplesResponses.rateZDuringPulse(all_ww)),1))...
    ,' SD'],'FontWeight','normal');

subplot(5,5,12)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
plotFill(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, rippleMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, rippleMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, rippleMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(rippleMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(rippleMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(rippleMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');
axis tight
xlabel('Ripple phase (deg)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% theta
subplot(5,5,13)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
plotFill(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, thetaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, thetaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, thetaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(thetaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(thetaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(thetaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');
axis tight
xlabel('Theta phase (deg)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% lgamma
subplot(5,5,14)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
plotFill(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, lgammaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, lgammaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, lgammaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(lgammaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(lgammaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(lgammaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');
axis tight
xlabel('Low gamma phase (deg)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% hgamma
subplot(5,5,15)
hold on
x_wave = 0:0.01:4*pi;
y_wave = cos(x_wave)*2;
plot(x_wave,y_wave,'--k');
plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_nw,:),'color',nw_color,'style','filled');
plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_ww,:),'color',ww_color,'style','filled');
plotFill(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(all_pyr,:),'color',pyr_color,'style','filled');

scatter(rand(length(find(all_pyr)),1)/2 + 12.6, hgammaMod.phasestats.r(all_pyr)*5,20,pyr_color,'filled');
scatter(rand(length(find(all_nw)),1)/2 + 13.1, hgammaMod.phasestats.r(all_nw)*5,20,nw_color,'filled');
scatter(rand(length(find(all_ww)),1)/2 + 13.6, hgammaMod.phasestats.r(all_ww)*5,20,ww_color,'filled');

scatter(hgammaMod.phasestats.m(all_pyr), rand(length(find(all_pyr)),1)/2 + 2.2,20,pyr_color,'filled');
scatter(hgammaMod.phasestats.m(all_nw), rand(length(find(all_nw)),1)/2 + 2.8,20,nw_color,'filled');
scatter(hgammaMod.phasestats.m(all_ww), rand(length(find(all_ww)),1)/2 + 3.4,20,ww_color,'filled');
axis tight
xlabel('High gamma phase (deg)'); ylabel('Rate (SD)');
set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});

% speed
subplot(5,5,16)
if ~isempty(speedCorr)
    if size(speedVals,1) ~= length(all_nw)
        speedVals = speedVals';
    end
    plotFill(log10(speedCorr.prc_vals),speedVals(all_nw,:),'color',nw_color,'style','filled');
    plotFill(log10(speedCorr.prc_vals),speedVals(all_ww,:),'color',ww_color,'style','filled');
    plotFill(log10(speedCorr.prc_vals),speedVals(all_pyr,:),'color',pyr_color,'style','filled');
    LogScale('x',10);

    scatter(rand(length(find(all_pyr)),1)/5 + 1.3, 1+speedCorr.speedScore(all_pyr)*2,20,pyr_color,'filled');
    scatter(rand(length(find(all_nw)),1)/5 +  1.4, 1+speedCorr.speedScore(all_nw)*2,20,nw_color,'filled');
    scatter(rand(length(find(all_ww)),1)/5 + 1.5, 1+speedCorr.speedScore(all_ww)*2,20,ww_color,'filled');
    axis tight
end
xlabel('cm/s'); ylabel('speed rate/avg rate');

subplot(5,5,17)
if ~isempty(behavior)
    hold on
    t_win = behavior.psth_reward.timestamps > -2 & behavior.psth_reward.timestamps < 2;
    plotFill(behavior.psth_reward.timestamps(t_win),behavior.psth_reward.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
    plotFill(behavior.psth_reward.timestamps(t_win),behavior.psth_reward.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
    plotFill(behavior.psth_reward.timestamps(t_win),behavior.psth_reward.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
    axis tight
    xlabel('Reward time (s)'); ylabel('Rate (SD)');

    subplot(5,5,21)
    hold on
    t_win = behavior.psth_intersection.timestamps > -2 & behavior.psth_intersection.timestamps < 2;
    plotFill(behavior.psth_intersection.timestamps(t_win),behavior.psth_intersection.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
    plotFill(behavior.psth_intersection.timestamps(t_win),behavior.psth_intersection.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
    plotFill(behavior.psth_intersection.timestamps(t_win),behavior.psth_intersection.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
    axis tight
    xlabel('Intersection time (s)'); ylabel('Rate (SD)');

    subplot(5,5,22)
    hold on
    t_win = behavior.psth_startPoint.timestamps > -2 & behavior.psth_startPoint.timestamps < 2;
    plotFill(behavior.psth_startPoint.timestamps(t_win),behavior.psth_startPoint.responsecurveZSmooth(all_nw,t_win),'color',nw_color,'style','filled');
    plotFill(behavior.psth_startPoint.timestamps(t_win),behavior.psth_startPoint.responsecurveZSmooth(all_ww,t_win),'color',ww_color,'style','filled');
    plotFill(behavior.psth_startPoint.timestamps(t_win),behavior.psth_startPoint.responsecurveZSmooth(all_pyr,t_win),'color',pyr_color,'style','filled');
    axis tight
    xlabel('Start point time (s)'); ylabel('Rate (SD)');
end

% spatial modulation
subplot(5,5,[18 23])
if ~isempty(spatialModulation)
    hold on
    yyaxis left
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
    yyaxis right
    set(gca,'YTick',[]);


    subplot(5,5,[19 24])
    hold on
    yyaxis left
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
    yyaxis right
    colormap jet
end

% black space
subplot(5,5,25)
axis off
title([{session.animal.geneticLine;[basenameFromBasepath(pwd)]}],'FontWeight','normal');

saveas(gcf,['SummaryFigures\Summary.png']);
    


cd(prevPath);
end