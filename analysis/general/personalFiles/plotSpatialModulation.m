function [] = plotSpatialModulation(varargin)
%        [] = plotSpatialModulation(varargin)
%
% Display spatial modulation summary per cell
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
addParameter(p,'gridAnalysis',false,@islogical);
addParameter(p,'tint',true,@islogical);
addParameter(p,'numSubplots',[],@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
UID = p.Results.UID;
saveFigure = p.Results.saveFigure;
checkUnits = p.Results.checkUnits;
use_deltaThetaEpochs = p.Results.use_deltaThetaEpochs;
excludePlot = p.Results.excludePlot;
gridAnalysis = p.Results.gridAnalysis;
tint = p.Results.tint;
numSubplots = p.Results.numSubplots;

% Color rules
color_R99 = [.9 .0 .0];
color_R95 = [.0 .9 .0];
color_ns = [.5 .5 .5];

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
% waveforms
all_waveforms = zscore(reshape([cell_metrics.waveforms.filt{:}],...
    [length(cell_metrics.waveforms.time{1}) length(cell_metrics.waveforms.filt)]));
pyr_color = [1 .7 .7];
nw_color = [.3 .3 1];
ww_color = [.7 1 1];

% acg
acg_time = [-50 : 0.5 : 50];
acg = cell_metrics.acg.narrow;
acg = acg./sum(acg);

% acg PeakTime
targetFile = dir('*.ACGPeak.cellinfo.mat'); 
acgPeak = importdata(targetFile.name);

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
[averageCCG] = getAverageCCG;

% ripples responses
targetFile = dir('*.ripples_psth.cellinfo.mat'); ripplesResponses = importdata(targetFile.name);
targetFile = dir('*.ripple_120-200.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
ripples = rippleMasterDetector;

% theta and gamma/s
targetFile = dir('*.theta_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.lgamma_20-60.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.hgamma_60-100.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

% Theta for REM and RUN
targetFile = dir('*.thetaRun_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
targetFile = dir('*.thetaREM_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);

% speed
speedCorr = getSpeedCorr('force',false);
if ~isempty(speedCorr)
    for ii = 1:size(speedCorr.speedVals,1)
        speedVals(ii,:) = mean(speedCorr.speedVals(ii,:,:),3)/cell_metrics.firingRate(ii);
    end
end
    
% speedVals = bsxfun(@rdivide,mean(speedCorr.speedVals,3),cell_metrics.firingRate(:));

% spatial modulation
if tint
    targetFile = dir('*.spatialModulation_tint.cellinfo.mat'); 
    if ~isempty(targetFile)
        load(targetFile.name);
    else 
        spatialModulation = [];
    end
else
    targetFile = dir('*.spatialModulation.cellinfo.mat'); 
    if ~isempty(targetFile)
        load(targetFile.name);
    else 
        spatialModulation = [];
    end
end

% spatial modulation 2 halves
if tint
    targetFile = dir('*.spatialModulation2Halves_tint.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        spatialModulation2Halves = [];
    end
else
    targetFile = dir('*.spatialModulation2Halves.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        spatialModulation2Halves = [];
    end
end
    
% behavioural events

targetFile = dir('*.behavior.cellinfo.mat');
if ~isempty(targetFile)
    load(targetFile.name);
else
    behavior = [];
end

% firing maps
if tint
    targetFile = dir('*.firingMapsAvg_tint.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        firingMaps = [];
    end
else
    targetFile = dir('*.firingMapsAvg.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        firingMaps = [];
    end
end


% firing maps 2 halves
if tint
    targetFile = dir('*.firingMapsAvg2Halves_tint.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        firingMaps2Halves = [];
    end
else
    targetFile = dir('*.firingMapsAvg2Halves.cellinfo.mat');
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        firingMaps2Halves = [];
    end
end


session = loadSession;
if length(firingMaps.rateMaps{1}) == 1
    numSubplots = 5;
elseif length(firingMaps.rateMaps{1}) == 2
    numSubplots = 6;
end

for ii = 1:length(UID)
    % putative cell type
    if strcmpi(cell_metrics.putativeCellType{ii},'Pyramidal Cell')
        color = pyr_color;
    elseif strcmpi(cell_metrics.putativeCellType{ii},'Narrow Interneuron')
        color = nw_color;
    elseif strcmpi(cell_metrics.putativeCellType{ii},'Wide Interneuron')
        color = ww_color;
    end
         
    % Plotting
    gcf = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Waveform
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on;
    plot(cell_metrics.waveforms.time{1},all_waveforms(:,ii),'color',color);
    axis tight; xlabel('ms'); ylabel('Waveform amp (SD)');
    title(['Cell: ', num2str(ii)]);
    
    % ACG
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    plot(acg_time,acg(:,ii),'color',color);
    axis tight; xlabel('ms'); ylabel('ACG (prob)');
    title(['Cell type: ', num2str(cell_metrics.putativeCellType{ii})]);
    
    % ACG Peak
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    plot(acgPeak.acg_time,acgPeak.acg_smoothed_norm(:,ii),'color',color);
    XTick = [-2 -1 0 1];
    set(gca,'XTick',XTick);
    XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
    set(gca,'XTickLabel',XTickLabels);
    ylabel('logACG (prob)'); xlabel('Time(s)');
    axis tight;
    
    % Cell position
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on;
    scatter(cell_metrics.general.chanCoords.x+100, cell_metrics.general.chanCoords.y,10,[.9 .9 .9],"filled");
    scatter(cell_metrics.trilat_x(ii)+100, cell_metrics.trilat_y(ii),30,color,'x','LineWidth',2);
    ax = axis;
    ylabel('depth (\mum)'); xlabel('\mum')
    xlim([min(cell_metrics.general.chanCoords.x)+50 max(cell_metrics.general.chanCoords.x)+150]);
    ylim([min(cell_metrics.general.chanCoords.y)-100 max(cell_metrics.general.chanCoords.y)+100]);
    title(['Region: ', num2str(cell_metrics.brainRegion{ii})]);
    
    % Firing Rate Stability
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    plot(spikemat.timestamps/60, spikemat.data(:,ii),'color',color);
    set(gca,'YScale','log');
    title(['Stability: ' num2str(round(cell_metrics.firingRateInstability(ii),2))],'FontWeight','normal');
    ylim([.1 50]); xlim([0 spikemat.timestamps(end)/60]);
    xlabel('Time (min)'); ylabel('Rate (Hz)');
    
    % Firing rate states
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    plot(1:5,states_rate(ii,:),'color',color);
    set(gca,'YScale','log');
    set(gca,'XTick',[1:5],'XTickLabel',{'All','Wake','Run','NREM','REM'},'XTickLabelRotation',45);
    xlim([.8 5.2]); ylabel('Rate (Hz)');
    
    % States index
%     h = get(gcf); nextPlot = length(h.Children)+1;
%     subplot(5,5,nextPlot)
%     hold on
%     scatter(rand(length(find(ii)),1)/5 + .9, run_quiet_index(ii),20,color,'filled');
%     scatter(rand(length(find(ii)),1)/5 + 1.9, rem_nrem_index(ii),20,color,'filled');
%     scatter(rand(length(find(ii)),1)/5 + 2.9, cv(ii),20,color,'filled');
% 
%     plot([.8 3.4],[0 0],'-k');
%     set(gca,'TickDir','out','XTick',[1:3],'XTickLabel',{'run-quiet','rem-nrem','cv'},'XTickLabelRotation',45);
%     ylabel('Index (SD)'); xlim([.8 3.4]);
    
    % Average CCG
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    plot(averageCCG.timestamps,averageCCG.ZmedianCCG(ii,:),'color',color);
    axis tight
    xlabel('Time (s)'); ylabel('Rate (SD)');
    
    % Ripples
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on;
    t_win = ripplesResponses.timestamps > -0.25 & ripplesResponses.timestamps < 0.25;
    plot(ripplesResponses.timestamps(t_win),ripplesResponses.responsecurveZSmooth(ii,t_win),'color',color);
    title([num2str(round(ripplesResponses.rateZDuringPulse(ii),2)) ,' SD'],'FontWeight','normal');
    plot(ripples.rippleStats.maps.timestamps, zscore(mean(ripples.rippleStats.maps.ripples_raw)) +20,'-k')
    axis tight
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    % ripple phase
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on
    x_wave = 0:0.01:4*pi;
    y_wave = cos(x_wave)*2;
    plot(x_wave,y_wave,'--k');
    plot(rippleMod.phasebins_wide_doubled,rippleMod.phasedistro_wide_smoothZ_doubled(ii,:),'color',color);
    axis tight
    xlabel('Ripple phase (rad)'); ylabel('Rate (SD)');
    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});
    
    % theta
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on
    x_wave = 0:0.01:4*pi;
    y_wave = cos(x_wave)*2;
    plot(x_wave,y_wave,'--k');
    plot(thetaMod.phasebins_wide_doubled,thetaMod.phasedistro_wide_smoothZ_doubled(ii,:),'color',color);
    axis tight
    xlabel('Theta phase (rad)'); ylabel('Rate (SD)');
    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});
    
    % lgamma
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on
    x_wave = 0:0.01:4*pi;
    y_wave = cos(x_wave)*2;
    plot(x_wave,y_wave,'--k');
    plot(lgammaMod.phasebins_wide_doubled,lgammaMod.phasedistro_wide_smoothZ_doubled(ii,:),'color',color);
    axis tight
    xlabel('Low gamma phase (rad)'); ylabel('Rate (SD)');
    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});
    
    % hgamma
    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
    hold on
    x_wave = 0:0.01:4*pi;
    y_wave = cos(x_wave)*2;
    plot(x_wave,y_wave,'--k');
    plot(hgammaMod.phasebins_wide_doubled,hgammaMod.phasedistro_wide_smoothZ_doubled(ii,:),'color',color);
    axis tight
    xlabel('High gamma phase (rad)'); ylabel('Rate (SD)');
    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0', '2\pi', '4\pi'});
    
    % speed
    if ~isempty(speedCorr)
        subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
        plot(log10(speedCorr.prc_vals),speedVals(ii,:),'color',color);
        LogScale('x',10);
        axis tight
        xlabel('cm/s'); ylabel('speed rate/avg rate');
    end
    
    % Behavior
    if ~isempty(behavior)
        if isfield(behavior,'psth_lReward') && isstruct(behavior.psth_lReward)
            subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
            hold on;
            t_win = behavior.psth_lReward.timestamps > -2 & behavior.psth_lReward.timestamps < 2;
            plot(behavior.psth_lReward.timestamps(t_win),behavior.psth_lReward.responsecurveZSmooth(ii,t_win),'color',color);
            axis tight
            xlabel('lReward time (s)'); ylabel('Rate (SD)');       
        end
        
        if isfield(behavior,'psth_rReward') && isstruct(behavior.psth_rReward)
            subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
            hold on;
            t_win = behavior.psth_rReward.timestamps > -2 & behavior.psth_rReward.timestamps < 2;
            plot(behavior.psth_rReward.timestamps(t_win),behavior.psth_rReward.responsecurveZSmooth(ii,t_win),'color',color);
            axis tight
            xlabel('rReward time (s)'); ylabel('Rate (SD)');       
        end
        
        if isfield(behavior,'psth_reward') && isstruct(behavior.psth_reward)
            subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
            hold on;
            t_win = behavior.psth_reward.timestamps > -2 & behavior.psth_reward.timestamps < 2;
            plot(behavior.psth_reward.timestamps(t_win),behavior.psth_reward.responsecurveZSmooth(ii,t_win),'color',color);
            axis tight
            xlabel('reward time (s)'); ylabel('Rate (SD)');       
        end
    end
    
    targetFile = dir([session.general.name,'.Behavior.mat']);
    if ~isempty(targetFile)
        load(targetFile.name);
    else
        behavior = [];
    end

    % Spatial Modulation rateMaps
    if ~isempty(firingMaps) && ~any(ismember(excludePlot,{lower('spatialModulation')}))
        for jj = 1:length(firingMaps.rateMaps{ii})
%                 if length(firingMaps.cmBin) == 2
%                     timestamps = 0:round(firingMaps.cmBin{jj},1):round(firingMaps.cmBin{jj},1)*(length(firingMaps.rateMaps{1}{jj})-1);
%                 else
%                     timestamps = 0:round(firingMaps.cmBin,1):round(firingMaps.cmBin,1)*(length(firingMaps.rateMaps{1}{jj})-1);
%                 end
                timestamps = 0:round(firingMaps.cmBin{jj},1):round(firingMaps.cmBin{jj},1)*(length(firingMaps.rateMaps{1}{jj})-1);

                if ~isempty(firingMaps2Halves) && ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
                    for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                        subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                        if size(firingMaps2Halves.rateMaps{ii}{jj}{kk},1) > 1
                            hold on;
                            plot(firingMaps2Halves.pos{jj}{kk}(:,2),firingMaps2Halves.pos{jj}{kk}(:,3),'color',[0.7 0.7 0.7]);
                            t = firingMaps2Halves.pos{jj}{kk}(:,1);
                            dt = diff(t); dt(end+1) = dt(end); dt(dt > firingMaps2Halves.params.maxGap) = firingMaps2Halves.params.maxGap;
                            n = CountInIntervals(spikes.times{ii},[t t+dt]);
                            scatter(firingMaps2Halves.pos{jj}{kk}(n > 0 ,2),firingMaps2Halves.pos{jj}{kk}(n > 0,3),1,'MarkerEdgeColor',[1 0 0], 'MarkerFaceColor',[0.9 0 0]);
                            axis ij;
                            axis square;
%                             xlim(round(behavior.avFrame{jj}.xSize)); ylim(round(behavior.avFrame{jj}.ySize));
                            if kk == 1
                                title(['Map ', num2str(jj) ,' 1st Half']);
                            elseif kk == 2
                                title(['Map ', num2str(jj) ,' 2nd Half']);
                            end
                        end
                            
                    end

                    if size(firingMaps.rateMaps{ii}{jj},1) > 1
                        subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                        hold on;
                        plot(behavior.maps{jj}(:,2),behavior.maps{jj}(:,3),'color',[0.7 0.7 0.7]);
                        t = behavior.maps{jj}(:,1);
                        dt = diff(t); dt(end+1) = dt(end); de(dt > firingMaps.params.maxGap) = firingMaps.params.maxGap;
                        n = CountInIntervals(spikes.times{ii},[t t+dt]);
                        scatter(behavior.maps{jj}(n > 0 ,2),behavior.maps{jj}(n > 0,3),1,'MarkerEdgeColor',[1 0 0], 'MarkerFaceColor',[0.9 0 0]);
                        axis ij;
                        axis square;
%                         xlim(round(behavior.avFrame{jj}.xSize)); ylim(round(behavior.avFrame{jj}.ySize));
                        title(['Map ', num2str(jj)]);
                        try
                            if ~isempty(spatialModulation2Halves)
                                stability_var = ['stability_map_',num2str(jj),'_twoHalves'];
                                stability = spatialModulation2Halves.(stability_var){ii};
                                xlabel('cm'); ylabel(['Stability: ',num2str(round(stability,3))]);
                            end
                        catch
                            
                        end
                    
                    end
                else
                    subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                    plot(timestamps, firingMaps.rateMaps{ii}{jj},'color',color);
                    xlim(timestamps([1 end]));
                    xlabel('cm'); ylabel(['Map ',num2str(jj),' (-3 to 3 SD)']);
                end
        end
    end
    
    % 2D Rate Maps with randomization
   if ~isempty(firingMaps) && ~any(ismember(excludePlot,{lower('spatialModulation')}))
        for jj = 1:length(firingMaps.rateMaps{ii})
            if size(firingMaps.rateMaps{ii}{jj},1) > 1 % 2D spatial map
                if ~isempty(firingMaps2Halves) && ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
                    
                    for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                        subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                        imagesc(firingMaps2Halves.rateMaps{ii}{jj}{kk});
                        colormap(jet(15))
                        axis ij;
                        axis square;
                        xlim([1 round(size(firingMaps2Halves.rateMaps{ii}{jj}{kk},1))]); ylim([1 round(size(firingMaps2Halves.rateMaps{ii}{jj}{kk},1))]);
                        if kk == 1
                            title(['Rate Map ', num2str(jj), ' 1st half']);
                        elseif kk == 2
                            title(['Rate Map ', num2str(jj), ' 2nd half']);
                        end
                    end
                    
                    ax = subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                    imagesc(firingMaps.rateMaps{ii}{jj});
                    colormap(jet(15))
                    axis ij;
                    axis square;
                    xlim([1 round(size(firingMaps.rateMaps{ii}{jj},1))]); 
                    ylim([1 round(size(firingMaps.rateMaps{ii}{jj},1))]);
                    
                    % Color Randomization
                    if tint
                        var_coherence = ['spatial_corr2_sc_r_map_',num2str(jj)];
                        coherence = spatialModulation.(var_coherence){ii};
                    else
                        var_coherence = ['spatial_corr_sc_r_map_',num2str(jj)];
                        coherence = spatialModulation.(var_coherence){ii};
                    end

                    var_bitsPerSpike = ['bitsPerSpike_map_',num2str(jj)];
                    bitsPerSpike = spatialModulation.(var_bitsPerSpike){ii};
                    
                    var_bitsPerSecond = ['bitsPerSec_map_',num2str(jj)];
                    bitsPerSecond = spatialModulation.(var_bitsPerSecond){ii};

                    shuffling_var = ['shuffling_map_',num2str(jj)];
                    
                    if tint
                        shuffling_coherenceR99 = spatialModulation.(shuffling_var){ii}.spatial_corr2_sc_r.R99;
                        shuffling_coherenceR95 = spatialModulation.(shuffling_var){ii}.spatial_corr2_sc_r.R95;
                    else
                        shuffling_coherenceR99 = spatialModulation.(shuffling_var){ii}.spatial_corr_sc_r.R99;
                        shuffling_coherenceR95 = spatialModulation.(shuffling_var){ii}.spatial_corr_sc_r.R95;
                    end
                    
                    shuffling_bitsPerSpikeR99 = spatialModulation.(shuffling_var){ii}.bitsPerSpike.R99;
                    shuffling_bitsPerSpikeR95 = spatialModulation.(shuffling_var){ii}.bitsPerSpike.R95;

                    shuffling_bitsPerSecondR99 = spatialModulation.(shuffling_var){ii}.bitsPerSec.R99;
                    shuffling_bitsPerSecondR95 = spatialModulation.(shuffling_var){ii}.bitsPerSec.R95;
                    
                    if coherence > shuffling_coherenceR99
                        color_coherence = color_R99;
                    elseif coherence > shuffling_coherenceR95
                        color_coherence = color_R95;
                    else
                        color_coherence = color_ns;
                    end

                    if bitsPerSpike > shuffling_bitsPerSpikeR99
                        color_bitsPerSpike = color_R99;
                    elseif bitsPerSpike > shuffling_bitsPerSpikeR95
                        color_bitsPerSpike = color_R95;
                    else
                        color_bitsPerSpike = color_ns;
                    end
                    
                    if bitsPerSecond > shuffling_bitsPerSecondR99
                        color_bitsPerSecond = color_R99;
                    elseif bitsPerSecond > shuffling_bitsPerSecondR95
                        color_bitsPerSecond = color_R95;
                    else
                        color_bitsPerSecond = color_ns;
                    end
                    colorbar;
                    yyaxis left
                    ylabel(ax,['Coherence: ', num2str(round(coherence,3))],'Color',color_coherence);
                    yyaxis right
                    set(gca,'ytick',[])
                    ylabel(ax,['BitsPerSecond: ', num2str(round(bitsPerSecond,3))],'Color',color_bitsPerSecond);
                    xlabel(ax,['BitsPerSpike: ' num2str(round(bitsPerSpike,3))],'color',color_bitsPerSpike);
                    firingFieldSize_var = ['firingFieldSize_map_',num2str(jj)];
                    title(['Max FR: ', num2str(spatialModulation.(firingFieldSize_var){ii}.maxFr)]); 
                end
            end
        end
   end
   
    % Spatial AutoCorrelation
    
    if ~isempty(spatialModulation)
        for jj = 1:length(firingMaps.rateMaps{ii})
            if size(firingMaps.rateMaps{ii}{jj},1) > 1 % 2D spatial map
                subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                spatialAutoCorr_var = ['spatialAutoCorr_map_',num2str(jj)];
                spatialAutoCorr = spatialModulation.(spatialAutoCorr_var){ii}.r;
                imagesc(spatialAutoCorr);
                colormap(jet(15))
                axis ij;
                axis square;
                xlim([1 round(size(spatialAutoCorr,1))]); ylim([1 round(size(spatialAutoCorr,2))]);
                title('Spatial AutoCorrelogram');
            end
        end
    end
    
    % Grid Analysis
    
    if ~isempty(spatialModulation) && gridAnalysis
        for jj = 1:length(firingMaps.rateMaps{ii})
            if size(firingMaps.rateMaps{ii}{jj},1) > 1 % 2D spatial map
                subplot(numSubplots,numSubplots,length(findobj(gcf,'type','axes'))+1)
                grid_var = ['grid_map_',num2str(jj)];
                gridMap = spatialModulation.(grid_var){ii};
                autoCorr = gridMap.autoCorr;
                if ~isnan(autoCorr.autoCorrMap)
                    plot(autoCorr.grad,autoCorr.autoCorrMap)
                    axis([0 180 -1 1]);
                    grid on;
                    axis square;
                    set(gca,'XTick',0:30:180,'XGrid','on','YGrid','on')
                    xlabel('Rotation [°]');
                    ylabel('autoCorrMap (r)');
                    if autoCorr.squareIndex > autoCorr.gridness
                        title(['Square Index: ', num2str(autoCorr.squareIndex)]);
                    elseif autoCorr.gridness > autoCorr.squareIndex
                        title(['Gridness Score: ' ,num2str(autoCorr.gridness)]);
                    end
                end
            end
        end
    end
    
    if saveFigure
        mkdir('spatialModulation');
        if tint
            saveas(gcf,['spatialModulation\Summary_Cell',num2str(ii),'_tint.png']);
        else
            saveas(gcf,['spatialModulation\Summary_Cell',num2str(ii),'.png']);
        end
    end
    close all;
end

cd(prevPath);
end




