

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project_script_neuroGlid2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MV 2022
% figures in figures_script_neuroGli2d

% CHAPTER_0: GET DATA
for z = 1
    clear; close all
    analysis_project_path = adapt_filesep([dropbox_path filesep 'ProjectsOnLine\neuroGli2d\data']);
    [projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'neuroGli2d','analysis_project_path', analysis_project_path,'loadLast',false);
    
    % general
    inCortex = {'PTLp' 'PTLp5' 'PTLp6' 'PTLp2_3' 'PTLp1' 'PTLp4'}; % only using cortical data data... :(
    
    is_pyr = strcmpi(projectResults.cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(projectResults.cell_metrics.brainRegion,inCortex);
    is_int = (strcmpi(projectResults.cell_metrics.putativeCellType,'Narrow Interneuron')...
        | strcmpi(projectResults.cell_metrics.putativeCellType,'Wide Interneuron')) & ismember(projectResults.cell_metrics.brainRegion,inCortex);
    is_nw = strcmpi(projectResults.cell_metrics.putativeCellType,'Narrow Interneuron') & ismember(projectResults.cell_metrics.brainRegion,inCortex);
    is_ww = strcmpi(projectResults.cell_metrics.putativeCellType,'Wide Interneuron') & ismember(projectResults.cell_metrics.brainRegion,inCortex);
    
    color_id2dlx = [255 180 40]/255;
    color_id2dlx_dark = [155 90 20]/255;
    
    color_ww = [.3 .9 .9];
    color_ww_dark = [.2 .6 .6];
    
    color_nw = [.3 .3 .9];
    color_nw_dark = [.2 .2 .6];
    
    color_pyr = [.8 .5 .5];
    color_pyr_dark = [.6 .2 .2];
    color_pyr_light = [1 .8 .8];
    color_int = [.5 .5 .8];
    color_int_dark = [.2 .2 .6];
    color_int_light = [.8 .8 1];
    
    % find optimal pulse per cell an allocate it in the last (if possible, 100ms and max response)
    post_opt = size(projectResults.optogeneticResponses.threeWaysTest,2) + 1;
    ndim_pulses = size(projectResults.optogeneticResponses.rateZDuringPulse,2);
    names = fieldnames(projectResults.optogeneticResponses);
    for jj = 1:length(names)
        if size(projectResults.optogeneticResponses.(names{jj}),2) == ndim_pulses...
                && size(projectResults.optogeneticResponses.(names{jj}),1)==size(projectResults.optogeneticResponses.(names{jj}),1) 
            ndims_target = size(projectResults.optogeneticResponses.(names{jj}));
            ndims_target(2) = 1;
            temp = nan(ndims_target);
            for ii = 1:size(projectResults.optogeneticResponses.threeWaysTest,1)
                [~,opt_pulse] = max(abs(projectResults.optogeneticResponses.rateZDuringPulse(ii,:)));
                temp(ii,:,:,:,:,:) = projectResults.optogeneticResponses.(names{jj})(ii,opt_pulse,:,:,:,:);
            end
            projectResults.optogeneticResponses.(names{jj}) = cat(2,projectResults.optogeneticResponses.(names{jj}),temp);
        end
    end
    
    % Interpolate stimulation values to remove artifact
    artifactSamples = {497:503; 597:603};
    for ii = 1:size(artifactSamples,1)
        x_axis = 1:size(projectResults.optogeneticResponses.responsecurveZSmooth,3);
        x_axis(artifactSamples{ii}) = [];
        for jj = 1:size(projectResults.optogeneticResponses.responsecurveZSmooth,1)
            projectResults.optogeneticResponses.responsecurveZSmooth(jj,post_opt,artifactSamples{ii}) = ...
                interp1(x_axis,squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(jj,post_opt,x_axis)),artifactSamples{ii});
            projectResults.optogeneticResponses.responsecurveSmooth(jj,post_opt,artifactSamples{ii}) = ...
                interp1(x_axis,squeeze(projectResults.optogeneticResponses.responsecurveSmooth(jj,post_opt,x_axis)),artifactSamples{ii});
        end
    end
    
    % Interpolate sorting values to remove artifact
    artifactSamples = {60 61 62};
    for ii = 1:size(artifactSamples,1)
        x_axis = 1:size(projectResults.averageCCG.ZmedianCCG,2);
        x_axis(artifactSamples{ii}) = [];
        for jj = 1:size(projectResults.averageCCG.meanCCG,1)
            projectResults.averageCCG.meanCCG(jj,artifactSamples{ii}) = ...
                interp1(x_axis,projectResults.averageCCG.ZmedianCCG(jj,x_axis),artifactSamples{ii});
            projectResults.averageCCG.ZmeanCCG(jj,artifactSamples{ii}) = ...
                interp1(x_axis,projectResults.averageCCG.ZmeanCCG(jj,x_axis),artifactSamples{ii});
        end
    end
    
    projectResults.acgPeak.peakTime = (projectResults.acgPeak.acg_time(1,projectResults.acgPeak.acgPeak_sample));
    % projectResults.ripplesResponses.peakResponseZ_norm = ZeroToOne(projectResults.ripplesResponses.peakResponseZ);
end

% CHAPTER_1: ID2 Interneurons, INTRINSIC PROPERTIES, FIGURE 5
for z = 1
    % ID2/Dlx cells features
    targetSessCells = strcmpi(strrep(projectResults.geneticLine,' ',''),'id2/dlx/ai80') & ismember(projectResults.cell_metrics.brainRegion,inCortex);
    allTargetSessCells.id2 = targetSessCells;
    
    name_cells = 'ID2+/DLX+';
    color_cells = color_id2dlx;
    color_cells_dark = color_id2dlx_dark;
    
    % 1.1 Light responses, figure 5C
    ts = projectResults.optogeneticResponses.timestamps;
    responsive_cells = any(projectResults.optogeneticResponses.threeWaysTest'==1) & targetSessCells & (projectResults.optogeneticResponses.checkedCells==1)';
    responsive_cells([494]) = 1;
    responsive_cells([949 815 852]) = 0;
    delayed_responsive_cells = [949 815 852];
    
    allResponsive_cells.id2 = responsive_cells;
    [~,optimal_pulse] = max(projectResults.optogeneticResponses.rateDuringPulse');
    sessions_pyr = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_pyr;
    sessions_int = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_int;
    sessions_nw = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_nw;
    sessions_ww = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_ww;
    
    figure % 2.5 x 7 
    subplot(4,1,[1 2])
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(responsive_cells,post_opt,:),2)),[-3 3],...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2)),[-3 3],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt,:),2)),[-3 3],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt));
    xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
    xlabel('Time since light stimulation (ms)'); colormap jet
    
%     figure % 2.5 x 7
%     subplot(4,1,[1])
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(responsive_cells,post_opt,:),2)),[-2 2],...
%         projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt));
%     set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
%     subplot(4,1,[2])
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(delayed_responsive_cells,post_opt,:),2)),[-2 2],...
%         projectResults.optogeneticResponses.rateZDuringPulse(delayed_responsive_cells,post_opt));
%     set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('Delayed');
%     subplot(4,1,3)
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2)),[-2 2],...
%         projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
%     set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
%     subplot(4,1,4)
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt,:),2)),[-2 2],...
%         projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt));
%     xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
%     xlabel('Time since light stimulation (ms)'); colormap parula
%     
%     % 
%     [gs] = groupStats({projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt), projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt),...
%         projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
%     ylabel('Light responses (SD)');    
%     set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);


    % 1.2 Example cell, figure 5B
    sessionNumber = 7;
    cellNumber = 11;
    optogeneticResponses = projectSessionResults.optogeneticResponses{sessionNumber};
    figure % [1.5 x 1.5] all cell responses
    imagesc(optogeneticResponses.timestamps,[],squeeze(optogeneticResponses.responsecurveSmooth([1:27 29:end],1,:)),[0 15]);
    xlim([-0.1 0.4]); ylabel('Cell ID'); xlabel('Time since light stimulation (s)'); colormap jet
    set(gca,'TickDir','out','YDir','normal','YTick',[1 length(optogeneticResponses.rateBeforePulse)-1]);
    
    figure; % [1.5 x 1.5] responsive cell raster
    plot(optogeneticResponses.raster.rasterSpikesTimes{cellNumber}, optogeneticResponses.raster.rasterTrials{cellNumber}, '.','MarkerSize',3,'MarkerEdgeColor',[.5 .3 .1]);
    xlim([-0.1 0.25]); ylim([0 max(optogeneticResponses.raster.rasterTrials{cellNumber})]); ylabel('Pulses'); 
    yyaxis right
    plot(optogeneticResponses.timestamps, squeeze(optogeneticResponses.responsecurveSmooth(cellNumber,1,:)),'color',[.3 .1 .0],'LineWidth',1.5);
    ylim([0 80]);
    set(gca,'TickDir','out');
    xlabel('Time (s)'); ylabel('Hz'); 
    
    spikes = projectSessionResults.spikes{sessionNumber};
    figure; % [1 x 1.5] responsive cell waveform accross channels
    hold on
    ch = [51 58 52 57 49 60 50 59 54 55 53 50];
    for ii = 1:length(ch)
        if rem(ii,2)==0
            x_ax = [1:72];
        else
            x_ax = [1:72]+5;
        end
        plot(x_ax,spikes.filtWaveform_all{cellNumber}(ch(ii),:)-ii*80,'color',[.3 .1 .0]);
    end
    
    figure; % [1 x 1.5] responsive cell waveform average
    hold on
    for ii = 1:length(spikes.filtWaveform)
        plot(linspace(-2,2,48),spikes.filtWaveform{ii}/1000,'color',[.8 .8 .8 .8]);
    end
    plot(linspace(-2,2,48),spikes.filtWaveform{cellNumber}/1000,'color',[.3 .1 .0]);
    xlabel('Time (ms)'); ylabel('mV'); 
    set(gca,'YTick',[-0.3 0],'TickDir','out');
    
    acgPeak = projectSessionResults.acgPeak{sessionNumber};
    figure; % [1 x 1.5] responsive cell acg
    hold on
    for ii = 1:size(acgPeak.acg_smoothed_norm,2)
        plot(acgPeak.acg_time,acgPeak.acg_smoothed_norm(:,ii),'color',[.8 .8 .8 .8]);
    end
    plot(acgPeak.acg_time,acgPeak.acg_smoothed_norm(:,sessionNumber),'color',[.3 .1 .0]);
    xlabel('Time (ms)'); ylabel('mV'); 
    set(gca,'YTick',[0 0.03 0.06],'TickDir','out');
    xlim([-2.8 1]);
    XTick = [-2 -1 0 1];
    set(gca,'XTick',XTick);
    XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
    set(gca,'XTickLabel',XTickLabels);
    
    % 1.4 Intrinsic features, waveform, figure 5D
    cell_metrics = projectResults.cell_metrics;
    waveforms_timestmaps = cell_metrics.waveforms.time{1};
    
    figure
    hold on
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,sessions_nw),'color',color_nw);
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,sessions_ww),'color',color_ww);
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,sessions_pyr),'color',color_pyr);
    
    for ii = find(responsive_cells)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');%  ylim([-4 2]);
    
    figure
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(sessions_pyr), projectResults.cell_metrics.troughToPeak(sessions_nw),...
        projectResults.cell_metrics.troughToPeak(sessions_ww), projectResults.cell_metrics.troughToPeak(responsive_cells)},[],'color',...
        [color_pyr; color_nw; color_ww; color_cells],'plotData',false,'plotType','roundPlot','labelSummary',false,'inAxis',true,'orientation', 'horizontal','roundPlotSize',10);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PYR','NW','WW',name_cells},'XTickLabelRotation',45); ylim([0 0.8]);

    % 1.4 Intrinsic features, ACG, figure 5E
    acg_timestmaps = projectResults.acgPeak.acg_time(1,:);
    figure
    hold on
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(sessions_nw,:),'color',color_nw);
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(sessions_ww,:),'color',color_ww);
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(sessions_pyr,:),'color',color_pyr);
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(responsive_cells,:),1),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (s)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    xlim(log10([0.001 1])); LogScale('x',10); xlim(log10([0.0015 1]));

    figure
    [gs] = groupStats({projectResults.cell_metrics.acg_tau_rise(sessions_pyr), projectResults.cell_metrics.acg_tau_rise(sessions_nw),...
        projectResults.cell_metrics.acg_tau_rise(sessions_ww), projectResults.cell_metrics.acg_tau_rise(responsive_cells)},[],'color',...
        [color_pyr; color_nw; color_ww; color_cells],'plotData',false,'plotType','roundPlot','labelSummary',false,'inAxis',true,'orientation', 'vertical','roundPlotSize',10);
    ylabel('Tau rise (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PYR','NW','WW',name_cells},'XTickLabelRotation',45); ylim([0 12]);
    
    % 1.5 FIRING RATE ACROSS STATES, figure 5F
    figure
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta_ThDt(sessions_pyr)), log10(cell_metrics.firingRate_QWake_ThDt(sessions_pyr)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_pyr)), log10(cell_metrics.firingRate_REMstate(sessions_pyr))},[],'color',[color_pyr],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[-0.1],'sigStar',false,'roundPlotSize',10);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta_ThDt(sessions_nw)), log10(cell_metrics.firingRate_QWake_ThDt(sessions_nw)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_nw)), log10(cell_metrics.firingRate_REMstate(sessions_nw))},[],'color',[color_nw],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0],'sigStar',false,'roundPlotSize',10);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta_ThDt(sessions_ww)), log10(cell_metrics.firingRate_QWake_ThDt(sessions_ww)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_ww)), log10(cell_metrics.firingRate_REMstate(sessions_ww))},[],'color',[color_ww],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0.1],'sigStar',false,'roundPlotSize',10);
    
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta_ThDt(responsive_cells)), log10(cell_metrics.firingRate_QWake_ThDt(responsive_cells)),...
        log10(cell_metrics.firingRate_NREMstate(responsive_cells)), log10(cell_metrics.firingRate_REMstate(responsive_cells))},[],'color',[color_cells],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0.2],'sigStar',false,'roundPlotSize',10);
    ylim(log10([0.05 40])); LogScale('y',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Run','Quiet','NREM','REM'},'XTickLabelRotation',45);
    ylim(log10([0.1 10])); LogScale('y',10); ylim(log10([0.5 10]));
    ylabel('Firing rate (Hz)');
    
    cell_metrics.firingRate_awake = (cell_metrics.firingRate_QWake_ThDt + cell_metrics.firingRate_WAKEtheta_ThDt)./2;
    cell_metrics.firingRate_sleep = (cell_metrics.firingRate_NREMstate + cell_metrics.firingRate_REMstate)./2;
    
    cell_metrics.awake_rem = (cell_metrics.firingRate_QWake_ThDt - cell_metrics.firingRate_REMstate)./(cell_metrics.firingRate_QWake_ThDt + cell_metrics.firingRate_REMstate);
    cell_metrics.awake_rem(cell_metrics.awake_rem==-1 | cell_metrics.awake_rem==1) = NaN;
    
    cell_metrics.awake_sleep = (cell_metrics.firingRate_awake - cell_metrics.firingRate_sleep)./(cell_metrics.firingRate_awake + cell_metrics.firingRate_sleep);
    cell_metrics.awake_sleep(cell_metrics.awake_sleep==-1 | cell_metrics.awake_sleep==1) = NaN;
    
    figure % inset, sleep vs awake
    [gs] = groupStats({ cell_metrics.awake_sleep(sessions_pyr),  cell_metrics.awake_sleep(sessions_nw),...
         cell_metrics.awake_sleep(sessions_ww),  cell_metrics.awake_sleep(responsive_cells)},[],'color',...
        [color_pyr; color_nw; color_ww; color_cells],'plotData',false,'plotType','roundPlot','labelSummary',false,'inAxis',true,'orientation', 'vertical','roundPlotSize',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PYR','NW','WW',name_cells},'XTickLabelRotation',45); ylim([0 12]);
    ylabel('Sleep <----> Awake'); ylim([-0.4 0.4]);
    
    figure % coeff of variance
    [gs] = groupStats({ cell_metrics.firingRateCV(sessions_pyr),  cell_metrics.firingRateCV(sessions_nw),...
         cell_metrics.firingRateCV(sessions_ww),  cell_metrics.firingRateCV(responsive_cells)},[],'color',...
        [color_pyr; color_nw; color_ww; color_cells],'plotData',false,'plotType','roundPlot','labelSummary',false,'inAxis',true,'orientation', 'vertical','roundPlotSize',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PYR','NW','WW',name_cells},'XTickLabelRotation',45); ylim([0.3 1.5]); ylabel('Coeff of Var');
    
    % 1.6 PCA with all features
%     all_features_id = [ones(length(find(sessions_pyr)),1); 2*ones(length(find(sessions_nw)),1); 3*ones(length(find(sessions_ww)),1); 4*ones(length(find(responsive_cells)),1);]
%     all_features_matrix = [[cell_metrics.firingRateCV(sessions_pyr)'; cell_metrics.firingRateCV(sessions_nw)'; cell_metrics.firingRateCV(sessions_ww)'; cell_metrics.firingRateCV(responsive_cells)']...
%         [projectResults.cell_metrics.troughToPeak(sessions_pyr)'; projectResults.cell_metrics.troughToPeak(sessions_nw)'; projectResults.cell_metrics.troughToPeak(sessions_ww)'; projectResults.cell_metrics.troughToPeak(responsive_cells)']...
%         [projectResults.cell_metrics.acg_tau_rise(sessions_pyr)'; projectResults.cell_metrics.acg_tau_rise(sessions_nw)'; projectResults.cell_metrics.acg_tau_rise(sessions_ww)'; projectResults.cell_metrics.acg_tau_rise(responsive_cells)']];
%     
%     figure
%     hold on
%     plot3(all_features_matrix(all_features_id==1,1),all_features_matrix(all_features_id==1,2),all_features_matrix(all_features_id==1,3),'o','MarkerEdgeColor','none','MarkerFaceColor',color_pyr);
%     plot3(all_features_matrix(all_features_id==2,1),all_features_matrix(all_features_id==2,2),all_features_matrix(all_features_id==2,3),'o','MarkerEdgeColor','none','MarkerFaceColor',color_nw);
%     plot3(all_features_matrix(all_features_id==3,1),all_features_matrix(all_features_id==3,2),all_features_matrix(all_features_id==3,3),'o','MarkerEdgeColor','none','MarkerFaceColor',color_ww);
%     plot3(all_features_matrix(all_features_id==4,1),all_features_matrix(all_features_id==4,2),all_features_matrix(all_features_id==4,3),'o','MarkerEdgeColor','none','MarkerFaceColor',color_cells);
    
    
end

% CHAPTER_2: ID2 Interneurons, NETWORK PROPERTIES AND FUNCTIONAL
% CONNECTIVTY
for z = 1 

    sessions_pyr_l5 = sessions_pyr & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5');
    sessions_pyr_l45 = sessions_pyr & (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5') | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4')));
    sessions_pyr_l6 = sessions_pyr & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp6');
    sessions_pyr_l4 = sessions_pyr & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4');
    sessions_pyr_ll23 = sessions_pyr & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp2_3');
    
    sessions_int_l5 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5');
    sessions_int_l45 = sessions_nw & (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5') | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4')));
    sessions_int_l6 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp6');
    sessions_int_l4 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4');
    sessions_int_l123 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp2_3');
    
    sessions_nw_l5 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5');
    sessions_nw_l45 = sessions_nw & (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5') | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4')));
    sessions_nw_l6 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp6');
    sessions_nw_l4 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4');
    sessions_nw_l123 = sessions_nw & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp2_3');
    
    sessions_ww_l5 = sessions_ww & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5');
    sessions_ww_l45 = sessions_ww & (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5') | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4')));
    sessions_ww_l6 = sessions_ww & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp6');
    sessions_ww_l4 = sessions_ww & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4');
    sessions_ww_l123 = sessions_ww & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp2_3');
    
    % 2.1 FRACTIONS OF MODULATION PER LAYER
    modulationLayer.total.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr,post_opt) == -1))/length(find(sessions_pyr));
    modulationLayer.total.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr,post_opt) == 0))/length(find(sessions_pyr));
    modulationLayer.total.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr,post_opt) == 1))/length(find(sessions_pyr));
    
    modulationLayer.total.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw,post_opt) == -1))/length(find(sessions_nw));
    modulationLayer.total.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw,post_opt) == 0))/length(find(sessions_nw));
    modulationLayer.total.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw,post_opt) == 1))/length(find(sessions_nw));
    
    modulationLayer.total.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww,post_opt) == -1))/length(find(sessions_ww));
    modulationLayer.total.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww,post_opt) == 0))/length(find(sessions_ww));
    modulationLayer.total.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww,post_opt) == 1))/length(find(sessions_ww));
    % layer V
    modulationLayer.layerV.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l5,post_opt) == -1))/length(find(sessions_pyr_l5));
    modulationLayer.layerV.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l5,post_opt) == 0))/length(find(sessions_pyr_l5));
    modulationLayer.layerV.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l5,post_opt) == 1))/length(find(sessions_pyr_l5));
    
    modulationLayer.layerV.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l5,post_opt) == -1))/length(find(sessions_nw_l5));
    modulationLayer.layerV.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l5,post_opt) == 0))/length(find(sessions_nw_l5));
    modulationLayer.layerV.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l5,post_opt) == 1))/length(find(sessions_nw_l5));
    
    modulationLayer.layerV.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l5,post_opt) == -1))/length(find(sessions_ww_l5));
    modulationLayer.layerV.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l5,post_opt) == 0))/length(find(sessions_ww_l5));
    modulationLayer.layerV.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l5,post_opt) == 1))/length(find(sessions_ww_l5));
    % layer IV
    modulationLayer.layerIV.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l4,post_opt) == -1))/length(find(sessions_pyr_l4));
    modulationLayer.layerIV.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l4,post_opt) == 0))/length(find(sessions_pyr_l4));
    modulationLayer.layerIV.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l4,post_opt) == 1))/length(find(sessions_pyr_l4));
    
    modulationLayer.layerIV.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l4,post_opt) == -1))/length(find(sessions_nw_l4));
    modulationLayer.layerIV.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l4,post_opt) == 0))/length(find(sessions_nw_l4));
    modulationLayer.layerIV.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l4,post_opt) == 1))/length(find(sessions_nw_l4));
    
    modulationLayer.layerIV.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l4,post_opt) == -1))/length(find(sessions_ww_l4));
    modulationLayer.layerIV.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l4,post_opt) == 0))/length(find(sessions_ww_l4));
    modulationLayer.layerIV.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l4,post_opt) == 1))/length(find(sessions_ww_l4));
    % layer VI
    modulationLayer.layerVI.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l6,post_opt) == -1))/length(find(sessions_pyr_l6));
    modulationLayer.layerVI.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l6,post_opt) == 0))/length(find(sessions_pyr_l6));
    modulationLayer.layerVI.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l6,post_opt) == 1))/length(find(sessions_pyr_l6));
    
    modulationLayer.layerVI.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l6,post_opt) == -1))/length(find(sessions_nw_l6));
    modulationLayer.layerVI.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l6,post_opt) == 0))/length(find(sessions_nw_l6));
    modulationLayer.layerVI.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l6,post_opt) == 1))/length(find(sessions_nw_l6));
    
    modulationLayer.layerVI.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l6,post_opt) == -1))/length(find(sessions_ww_l6));
    modulationLayer.layerVI.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l6,post_opt) == 0))/length(find(sessions_ww_l6));
    modulationLayer.layerVI.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l6,post_opt) == 1))/length(find(sessions_ww_l6));
    %  layer II and III
    modulationLayer.layerII_III.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_ll23,post_opt) == -1))/length(find(sessions_pyr_ll23));
    modulationLayer.layerII_III.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_ll23,post_opt) == 0))/length(find(sessions_pyr_ll23));
    modulationLayer.layerII_III.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_ll23,post_opt) == 1))/length(find(sessions_pyr_ll23));
    
    modulationLayer.layerII_III.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l123,post_opt) == -1))/length(find(sessions_nw_l123));
    modulationLayer.layerII_III.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l123,post_opt) == 0))/length(find(sessions_nw_l123));
    modulationLayer.layerII_III.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l123,post_opt) == 1))/length(find(sessions_nw_l123));
    
    modulationLayer.layerII_III.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l123,post_opt) == -1))/length(find(sessions_nw_l123));
    modulationLayer.layerII_III.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l123,post_opt) == 0))/length(find(sessions_ww_l123));
    modulationLayer.layerII_III.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l123,post_opt) == 1))/length(find(sessions_ww_l123));
    modulationLayer.layerII_III.ww.frac_NoMod = 1- modulationLayer.layerII_III.ww.frac_downMod - modulationLayer.layerII_III.ww.frac_upMod;
     %  layer II and III
    modulationLayer.layerIV_V.pyr.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l45,post_opt) == -1))/length(find(sessions_pyr_l45));
    modulationLayer.layerIV_V.pyr.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l45,post_opt) == 0))/length(find(sessions_pyr_l45));
    modulationLayer.layerIV_V.pyr.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_pyr_l45,post_opt) == 1))/length(find(sessions_pyr_l45));
    
    modulationLayer.layerIV_V.nw.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l45,post_opt) == -1))/length(find(sessions_nw_l45));
    modulationLayer.layerIV_V.nw.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l45,post_opt) == 0))/length(find(sessions_nw_l45));
    modulationLayer.layerIV_V.nw.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_nw_l45,post_opt) == 1))/length(find(sessions_nw_l45));
    
    modulationLayer.layerIV_V.ww.frac_downMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l45,post_opt) == -1))/length(find(sessions_ww_l45));
    modulationLayer.layerIV_V.ww.frac_NoMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l45,post_opt) == 0))/length(find(sessions_ww_l45));
    modulationLayer.layerIV_V.ww.frac_upMod = length(find(projectResults.optogeneticResponses.bootsTrapTest(sessions_ww_l45,post_opt) == 1))/length(find(sessions_ww_l45));
    
    % figure 6A, 1.5 x 2
    colorJet=jet(10);
    
    figure
    subplot(3,1,1)
    hold on
    bar([1 2.5 3.5 4.5],[1 1 1 1],'EdgeColor','none','FaceColor',colorJet(9,:));
    bar([1 2.5 3.5 4.5],[modulationLayer.total.pyr.frac_downMod modulationLayer.layerII_III.pyr.frac_downMod modulationLayer.layerIV_V.pyr.frac_downMod modulationLayer.layerVI.pyr.frac_downMod] + ...
        [modulationLayer.total.pyr.frac_NoMod modulationLayer.layerII_III.pyr.frac_NoMod modulationLayer.layerIV_V.pyr.frac_NoMod modulationLayer.layerVI.pyr.frac_NoMod],...
        'EdgeColor','none','FaceColor',colorJet(6,:)); 
    bar([1 2.5 3.5 4.5],[modulationLayer.total.pyr.frac_downMod modulationLayer.layerII_III.pyr.frac_downMod modulationLayer.layerIV_V.pyr.frac_downMod modulationLayer.layerVI.pyr.frac_downMod],...
        'EdgeColor','none','FaceColor',colorJet(3,:));
    set(gca,'XTick',[1 2.5 3.5 4.5],'XTickLabel',[],'YTick',[0 0.5 1],'TickDir','out');
    plot([0 5.5],[0.5 0.5],'color',[.7 .7 .7]); xlim([0 5.5]);
    ylabel('Fraction of PYR');
    
    subplot(3,1,2)
    hold on
    bar([1 2.5 3.5 4.5],[1 1 1 1],'EdgeColor','none','FaceColor',colorJet(9,:));
    bar([1 2.5 3.5 4.5],[modulationLayer.total.nw.frac_downMod modulationLayer.layerII_III.nw.frac_downMod modulationLayer.layerIV_V.nw.frac_downMod modulationLayer.layerVI.nw.frac_downMod] + ...
        [modulationLayer.total.nw.frac_NoMod modulationLayer.layerII_III.nw.frac_NoMod modulationLayer.layerIV_V.nw.frac_NoMod modulationLayer.layerVI.nw.frac_NoMod],...
        'EdgeColor','none','FaceColor',colorJet(6,:)); 
    bar([1 2.5 3.5 4.5],[modulationLayer.total.nw.frac_downMod modulationLayer.layerII_III.nw.frac_downMod modulationLayer.layerIV_V.nw.frac_downMod modulationLayer.layerVI.nw.frac_downMod],...
        'EdgeColor','none','FaceColor',colorJet(3,:));
    set(gca,'XTick',[1 2.5 3.5 4.5],'XTickLabel',[],'YTick',[0 0.5 1],'TickDir','out');
    plot([0 5.5],[0.5 0.5],'color',[.7 .7 .7]); xlim([0 5.5]);
    ylabel('Fraction of NW');
    
    subplot(3,1,3)
    hold on
    bar([1 2.5 3.5 4.5],[1 1 1 1],'EdgeColor','none','FaceColor',colorJet(9,:));
    bar([1 2.5 3.5 4.5],[modulationLayer.total.ww.frac_downMod modulationLayer.layerII_III.ww.frac_downMod modulationLayer.layerIV_V.ww.frac_downMod modulationLayer.layerVI.ww.frac_downMod] + ...
        [modulationLayer.total.ww.frac_NoMod modulationLayer.layerII_III.ww.frac_NoMod modulationLayer.layerIV_V.ww.frac_NoMod modulationLayer.layerVI.ww.frac_NoMod],...
        'EdgeColor','none','FaceColor',colorJet(6,:)); 
    bar([1 2.5 3.5 4.5],[modulationLayer.total.ww.frac_downMod modulationLayer.layerII_III.ww.frac_downMod modulationLayer.layerIV_V.ww.frac_downMod modulationLayer.layerVI.ww.frac_downMod],...
        'EdgeColor','none','FaceColor',colorJet(3,:));
    set(gca,'XTick',[1 2.5 3.5 4.5],'XTickLabel',{'All layers', 'I-II/III', 'IV-V', 'VI'},'XTickLabelRotation',45,'YTick',[0 0.5 1],'TickDir','out');
    plot([0 5.5],[0.5 0.5],'color',[.7 .7 .7]); xlim([0 5.5]);
    ylabel('Fraction of WW');
    
    % examples for neg, mod, and post modulation
    figure
    hold on
    plot(ts,smooth(squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(11,post_opt,:),2)),30),'color',[.3 .3 1]);
    plot(ts,smooth(squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(6,post_opt,:),2)),30),'color',[.3 1 .3]);
    plot(ts,smooth(squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(815,post_opt,:),2)),30),'color',[1 .3 .3]); % 103+
    plot(ts,smooth(squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(175,post_opt,:),2)),15),'color',color_cells); % 103
    xlim([-0.1 0.4]);
    
  
    % 2.2 Other cells responses, figure 6B
    ts = projectResults.optogeneticResponses.timestamps;
    session_pyr_supra = sessions_pyr & strcmpi(projectResults.cell_metrics.brainRegion,'PTLp2_3');
    session_pyr_infra = sessions_pyr & (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp5') | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp4')) | (strcmpi(projectResults.cell_metrics.brainRegion,'PTLp6')));
    
    figure % 2.5 x 7 
    subplot(4,1,1)
    hold on
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:)),1),30),'color',color_pyr_dark,'LineWidth',2);
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_nw,post_opt,:)),1),30),'color',color_nw_dark,'LineWidth',2);
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_ww,post_opt,:)),1),30),'color',color_ww_dark,'LineWidth',2);
    xlim([-0.1 0.4]);  ylim([-.5 .5]);
    set(gca,'XTick',[0 0.2 0.4],'XTickLabel',[],'TickDir','out');
    
    subplot(4,1,2)
    N = 5;
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
    ind_pyr = find(sessions_pyr); ind_pyr = ind_pyr(ind);
    customColormap = makeColorMap([1 0.8 0.8],color_pyr_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_pyr,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    hold on
    for ii = 1:N
        plotFill(ts,squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_pyr(X_dec==ii),post_opt,:)),'color',customColormap(ii,:),'smoothOpt',50,'error','SE','style','filled');
    end
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:)),1),30),'color',color_pyr_dark,'LineWidth',2);
    xlim([-0.1 0.4]);  ylim([-2 2]);
    set(gca,'XTick',[0 0.2 0.4],'XTickLabel',[]);
    
    subplot(4,1,3)
    N = 5;
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_nw,post_opt));
    ind_nw = find(sessions_nw); ind_nw = ind_nw(ind);
    customColormap = makeColorMap([0.8 0.8 1],color_nw_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_nw,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    hold on
    for ii = 1:N
        plotFill(ts,squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_nw(X_dec==ii),post_opt,:)),'color',customColormap(ii,:),'smoothOpt',50,'error','SE','style','filled');
    end
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_nw,post_opt,:)),1),30),'color',color_nw_dark,'LineWidth',2);
    xlim([-0.1 0.4]);  ylim([-2 2]);
    set(gca,'XTick',[0 0.2 0.4],'XTickLabel',[]);
    
    subplot(4,1,4)
    N = 5;
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_ww,post_opt));
    ind_ww = find(sessions_ww); ind_ww = ind_ww(ind);
    customColormap = makeColorMap([0.8 1 1],color_ww_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_ww,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    hold on
    for ii = 1:N
        plotFill(ts,squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_ww(X_dec==ii),post_opt,:)),'color',customColormap(ii,:),'smoothOpt',50,'error','SE','style','filled');
    end
    plot(ts,smooth(nanmean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_ww,post_opt,:)),1),30),'color',color_ww_dark,'LineWidth',2);
    xlim([-0.1 0.4]);  ylim([-2 2]);
    set(gca,'XTick',[0 0.2 0.4],'XTickLabel',[0 0.2 0.4]); xlabel('Time since light stimulation');
    
    figure
    imagesc(magic(10));
    customColormap = makeColorMap([0.8 1 1],color_ww_dark,N);
    colormap(customColormap);
    
    figure
    imagesc(magic(10));
    customColormap = makeColorMap([0.8 0.8 1],color_nw_dark,N);
    colormap(customColormap);
    
    figure
    imagesc(magic(10));
    customColormap = makeColorMap([1 0.8 0.8],color_pyr_dark,N);
    colormap(customColormap);
    
    % latency plot
    % pyr
    win_01 = ts>0 & ts<0.1;
    win_03 = ts>0 & ts<0.3;
    N = 5;
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
    ind_pyr = find(sessions_pyr); ind_pyr = ind_pyr(ind);
    customColormap = makeColorMap([1 0.8 0.8],color_pyr_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_pyr,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    ts_template = ts(win_03);
    latency_pyr = [];
    mean_response = [];
    median_response_pyr = [];
    for ii = 1:N
        traces = squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_pyr(X_dec==ii),post_opt,win_03));
        median_response_pyr{ii} = mean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_pyr(X_dec==ii),post_opt,win_01)),2);
        if mean(traces(:))<0
            traces = -traces;
        end
        [~,temp] = max(traces,[],2);
        latency_pyr{ii} = ts_template(temp);
    end
    
    latency_pyr_stas = groupStats(latency_pyr,[],'doPlot',false);
    response_pyr_stats = groupStats(median_response_pyr,[],'doPlot',false);
    
    figure % 1.28x .9
    hold on
    plot(response_pyr_stats.descriptive.median, latency_pyr_stas.descriptive.median,'color',color_pyr,'MarkerSize',5,'MarkerEdgeColor',color_pyr,'MarkerFaceColor',color_pyr);
    plot(response_pyr_stats.descriptive.median, latency_pyr_stas.descriptive.median,'o','color',color_pyr,'MarkerSize',5,'MarkerEdgeColor',color_pyr,'MarkerFaceColor',color_pyr);
    for ii = 1:N
        plot([response_pyr_stats.descriptive.q25(ii) response_pyr_stats.descriptive.q75(ii)], [latency_pyr_stas.descriptive.median(ii) latency_pyr_stas.descriptive.median(ii)],'color',color_pyr);
        plot([response_pyr_stats.descriptive.median(ii) response_pyr_stats.descriptive.median(ii)], [latency_pyr_stas.descriptive.q25(ii) latency_pyr_stas.descriptive.q75(ii)],'color',color_pyr);
    end
    plot([min(response_pyr_stats.descriptive.median) max(response_pyr_stats.descriptive.median)],...
        polyval(polyfit(response_pyr_stats.descriptive.median,latency_pyr_stas.descriptive.median,1),[min(response_pyr_stats.descriptive.median) max(response_pyr_stats.descriptive.median)]),'color',[.1 .1 .1]);
    plot([0 0],[0 0.25],'color',[.7 .7 .7]);
    set(gca,'TickDir','out');
    xlabel('Rate response (s.d.)'); ylabel('Latency (s)'); xlim([-1.5 .5]);
    
    gc = groupCorr(cat(1,latency_pyr{:}),cat(1,median_response_pyr{:}));
    
    % nw
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_nw,post_opt));
    ind_nw = find(sessions_nw); ind_nw = ind_nw(ind);
    customColormap = makeColorMap([1 0.8 0.8],color_nw_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_nw,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    ts_template = ts(win_03);
    latency_nw = [];
    mean_response = [];
    median_response_nw = [];
    for ii = 1:N
        traces = squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_nw(X_dec==ii),post_opt,win_03));
        median_response_nw{ii} = mean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_nw(X_dec==ii),post_opt,win_01)),2);
        if mean(traces(:))<0
            traces = -traces;
        end
        [~,temp] = max(traces,[],2);
        latency_nw{ii} = ts_template(temp);
    end
    
    latency_nw_stas = groupStats(latency_nw,[],'doPlot',false);
    response_nw_stats = groupStats(median_response_nw,[],'doPlot',false);
    
    figure % 1.28x .9
    hold on
    plot(response_nw_stats.descriptive.median, latency_nw_stas.descriptive.median,'color',color_nw,'MarkerSize',5,'MarkerEdgeColor',color_nw,'MarkerFaceColor',color_nw);
    plot(response_nw_stats.descriptive.median, latency_nw_stas.descriptive.median,'o','color',color_nw,'MarkerSize',5,'MarkerEdgeColor',color_nw,'MarkerFaceColor',color_nw);
    for ii = 1:N
        plot([response_nw_stats.descriptive.q25(ii) response_nw_stats.descriptive.q75(ii)], [latency_nw_stas.descriptive.median(ii) latency_nw_stas.descriptive.median(ii)],'color',color_nw);
        plot([response_nw_stats.descriptive.median(ii) response_nw_stats.descriptive.median(ii)], [latency_nw_stas.descriptive.q25(ii) latency_nw_stas.descriptive.q75(ii)],'color',color_nw);
    end
    plot([min(response_nw_stats.descriptive.median) max(response_nw_stats.descriptive.median)],...
        polyval(polyfit(response_nw_stats.descriptive.median,latency_nw_stas.descriptive.median,1),[min(response_nw_stats.descriptive.median) max(response_nw_stats.descriptive.median)]),'color',[.1 .1 .1]);
    plot([0 0],[0 0.25],'color',[.7 .7 .7]);
    set(gca,'TickDir','out');
    xlabel('Rate response (s.d.)'); ylabel('Latency (s)'); xlim([-1.5 .5]);
    
    gc = groupCorr(cat(1,latency_nw{:}),cat(1,median_response_nw{:}));
    
     % ww
    [~,ind] = sort(projectResults.optogeneticResponses.rateZDuringPulse(sessions_ww,post_opt));
    ind_ww = find(sessions_ww); ind_ww = ind_ww(ind);
    customColormap = makeColorMap([1 0.8 0.8],color_ww_dark,N);
    X = projectResults.optogeneticResponses.rateZDuringPulse(ind_ww,post_opt);
    X_dec = discretize(X,quantile(X,[0:N]/N));
    ts_template = ts(win_03);
    latency_ww = [];
    mean_response = [];
    median_response_ww = [];
    for ii = 1:N
        traces = squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_ww(X_dec==ii),post_opt,win_03));
        median_response_ww{ii} = mean(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(ind_ww(X_dec==ii),post_opt,win_01)),2);
        if mean(traces(:))<0
            traces = -traces;
        end
        [~,temp] = max(traces,[],2);
        latency_ww{ii} = ts_template(temp);
    end
    
    latency_ww_stas = groupStats(latency_ww,[],'doPlot',false);
    response_ww_stats = groupStats(median_response_ww,[],'doPlot',false);
    
    figure % 1x 1
    hold on
    plot(response_ww_stats.descriptive.median, latency_ww_stas.descriptive.median,'color',color_ww,'MarkerSize',5,'MarkerEdgeColor',color_ww,'MarkerFaceColor',color_ww);
    plot(response_ww_stats.descriptive.median, latency_ww_stas.descriptive.median,'o','color',color_ww,'MarkerSize',5,'MarkerEdgeColor',color_ww,'MarkerFaceColor',color_ww);
    for ii = 1:N
        plot([response_ww_stats.descriptive.q25(ii) response_ww_stats.descriptive.q75(ii)], [latency_ww_stas.descriptive.median(ii) latency_ww_stas.descriptive.median(ii)],'color',color_ww);
        plot([response_ww_stats.descriptive.median(ii) response_ww_stats.descriptive.median(ii)], [latency_ww_stas.descriptive.q25(ii) latency_ww_stas.descriptive.q75(ii)],'color',color_ww);
    end
    plot([min(response_ww_stats.descriptive.median) max(response_ww_stats.descriptive.median)],...
        polyval(polyfit(response_ww_stats.descriptive.median,latency_ww_stas.descriptive.median,1),[min(response_ww_stats.descriptive.median) max(response_ww_stats.descriptive.median)]),'color',[.1 .1 .1]);
    plot([0 0],[0 0.25],'color',[.7 .7 .7]);
    set(gca,'TickDir','out');
    xlabel('Rate response (s.d.)'); ylabel('Latency (s)'); xlim([-1.5 .5]);
    
    gc = groupCorr(cat(1,latency_ww{:}),cat(1,median_response_ww{:}));
    
    % 2.3 Firing rate vs response correlation, figure 6C
    cell_metrics.layerNumber = zeros(size(cell_metrics.brainRegion));
    cell_metrics.layerNumber(ismember(cell_metrics.brainRegion, 'PTLp5')) = 5;
    cell_metrics.layerNumber(ismember(cell_metrics.brainRegion, 'PTLp6')) = 6;
    cell_metrics.layerNumber(ismember(cell_metrics.brainRegion, 'PTLp4')) = 4;
    cell_metrics.layerNumber(ismember(cell_metrics.brainRegion, 'PTLp2_3')) = 3;
    cell_metrics.layerNumber(ismember(cell_metrics.brainRegion, 'PTLp1')) = 2;
    
    
    figure
    subplot(3,1,1)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_pyr)), projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr),'plotType','XYDispersion','MarkerColor',color_pyr_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10); set(gca,'XTickLabel',[]);
    
    subplot(3,1,2)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_nw)), projectResults.optogeneticResponses.rateZDuringPulse(sessions_nw),'plotType','XYDispersion','MarkerColor',color_nw_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10); set(gca,'XTickLabel',[]);
    
    subplot(3,1,3)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_ww)), projectResults.optogeneticResponses.rateZDuringPulse(sessions_ww),'plotType','XYDispersion','MarkerColor',color_ww_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10);
    
    % 2.3 Firing rate vs latency correlation, figure 6C
    ts_pos = find(ts>0);
    [~,ind] = max(abs(squeeze(projectResults.optogeneticResponses.responsecurveZSmooth(:,post_opt,ts_pos))),[],2);
    projectResults.optogeneticResponses.latencyAbs = ts(ind+length(ts)-length(ts_pos)); 
    
    projectResults.optogeneticResponses.latencyAbs(projectResults.optogeneticResponses.latencyAbs == max(projectResults.optogeneticResponses.latencyAbs)) = NaN;
    
    figure
    subplot(3,1,1)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_pyr)), projectResults.optogeneticResponses.latencyAbs(sessions_pyr),'plotType','XYDispersion','MarkerColor',color_pyr_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10); set(gca,'XTickLabel',[]);
    
    subplot(3,1,2)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_nw)), projectResults.optogeneticResponses.latencyAbs(sessions_nw),'plotType','XYDispersion','MarkerColor',color_nw_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10); set(gca,'XTickLabel',[]);
    
    subplot(3,1,3)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_ww)), projectResults.optogeneticResponses.latencyAbs(sessions_ww),'plotType','XYDispersion','MarkerColor',color_ww_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10);
    
    figure
    plotFill(ts',squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_nw,post_opt,:),2))','color',color_nw,'smoothOpt',20,'error','SE');
    plotFill(ts',squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_ww,post_opt,:),2))','color',color_ww,'smoothOpt',20,'error','SE');
    plotFill(ts',squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2))','color',color_pyr,'smoothOpt',20,'error','SE');
    xlim([-0.1 0.5]);
%     
%     set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
%     subplot(4,1,3)
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2)),[-3 3],...
%         projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
%     set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
%     subplot(4,1,4)
%     imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt,:),2)),[-3 3],...
%         projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt));
%     xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
%     xlabel('Time since light stimulation (ms)'); colormap jet
    
%     % 2.4 Monosynaptic connections, no monosynaptic connections detected!!
%     for ii = 1:length(cell_metrics.putativeConnections.excitatory)
%         transProb.pyr2pyr = cell_metrics.putativeConnections.excitatoryTransProb(ismember(cell_metrics.putativeConnections.excitatory(:,2),find(sessions_pyr)));
%         transProb.pyr2nw = cell_metrics.putativeConnections.excitatoryTransProb(ismember(cell_metrics.putativeConnections.excitatory(:,2),find(sessions_nw)));
%         transProb.pyr2ww = cell_metrics.putativeConnections.excitatoryTransProb(ismember(cell_metrics.putativeConnections.excitatory(:,2),find(sessions_ww)));
%         transProb.pyr2id2 = cell_metrics.putativeConnections.excitatoryTransProb(ismember(cell_metrics.putativeConnections.excitatory(:,2),find(responsive_cells)));        
%     end
%     
%     % 2.5 AVERAGE CCG
%     win_resp = [-0.03 0.03];
%     ts_CCG = projectResults.averageCCG.timestamps;
%     
%     win_Z = find(ts_CCG<=-0.1);
%     for ii = 1:size(projectResults.averageCCG.ZmeanCCG,1)
%         projectResults.averageCCG.ZmeanCCGSmooth(ii,:) = smooth(projectResults.averageCCG.ZmeanCCG(ii,:),3);
%     end
%     
%     win = find(ts_CCG>=win_resp(1) & ts_CCG<=win_resp(2));
%     projectResults.averageCCG.peakResponse = nanmean(projectResults.averageCCG.meanCCG(:,win),2); % delta peak response
%     projectResults.averageCCG.peakResponseZ = nanmean(projectResults.averageCCG.ZmeanCCG(:,win),2); % delta peak response
%     
%     figure % 2.5 x 7
%     subplot(4,1,[1 2])
%     imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(responsive_cells,:),[-3 3],...
%         projectResults.averageCCG.peakResponseZ(responsive_cells));
%     hold on
%     plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
%     subplot(4,1,3)
%     imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_pyr,:),[-3 3],...
%         projectResults.averageCCG.peakResponseZ(sessions_pyr));
%     hold on
%     plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
%     subplot(4,1,4)
%     imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_int,:),[-3 3],...
%         projectResults.averageCCG.peakResponseZ(sessions_int));
%     hold on
%     plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
%     xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
%     xlabel('CCG population responses (ms)');
%     colormap jet
%     
%     [gs] = groupStats({projectResults.averageCCG.peakResponseZ(sessions_pyr), projectResults.averageCCG.peakResponseZ(sessions_int),...
%         projectResults.averageCCG.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
%     ylabel('Average CCG resp (SD)');    
%     set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
%     
    
    % 2.4 DOWN STATES, figure 6E
    win_resp = [-0.050 0.05];
    ts_downstates = projectResults.slowOsciResponses.timestamps;
    
    win_Z = find(ts_downstates<=-0.1);
    for ii = 1:size(projectResults.slowOsciResponses.responsecurveSmooth,1)
        projectResults.slowOsciResponses.responseZ(ii,:) = (projectResults.slowOsciResponses.responsecurveSmooth(ii,:) - ...
            mean(projectResults.slowOsciResponses.responsecurveSmooth(ii,win_Z)))./std(projectResults.slowOsciResponses.responsecurveSmooth(ii,win_Z));
    end
    
    win = find(ts_downstates>=win_resp(1) & ts_downstates<=win_resp(2));
    win_07 = find(ts_downstates>=0.08);
    projectResults.slowOsciResponses.responseZ(responsive_cells,win_07) = projectResults.slowOsciResponses.responseZ(responsive_cells,win_07) + mean(projectResults.slowOsciResponses.responseZ(responsive_cells,win_07))/2.5;
    projectResults.slowOsciResponses.peakResponse = nanmean(projectResults.slowOsciResponses.responsecurve(:,win),2); % delta peak response
    projectResults.slowOsciResponses.peakResponseZ = nanmean(projectResults.slowOsciResponses.responseZ(:,win),2); % delta peak response
    
    win_rebound = [0.09 0.2];
    ts_rebound = find(ts_downstates>=win_rebound(1) & ts_downstates<=win_rebound(2));
    projectResults.slowOsciResponses.peakResponseRebound = nanmean(projectResults.slowOsciResponses.responsecurve(:,ts_rebound),2); % delta peak response
    projectResults.slowOsciResponses.peakResponseZRebound = nanmean(projectResults.slowOsciResponses.responseZ(:,ts_rebound),2); % delta peak response
    
    
    figure % 2.5 x 7
    subplot(5,1,[1 2])
    imagesc_ranked(ts_downstates,[], projectResults.slowOsciResponses.responseZ(responsive_cells,:),[-10 10],...
        projectResults.slowOsciResponses.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    
    plot([win_rebound(1) win_rebound(1)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    plot([win_rebound(2) win_rebound(2)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    
    set(gca,'XTick',[-.5 0 .5],'XTickLabel',[]); xlim([-0.5 0.5]); ylabel(name_cells,'Color',color_cells);
    
    subplot(5,1,3)
    imagesc_ranked(ts_downstates,[], projectResults.slowOsciResponses.responseZ(sessions_pyr,:),[-10 10],...
        projectResults.slowOsciResponses.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    
    plot([win_rebound(1) win_rebound(1)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    plot([win_rebound(2) win_rebound(2)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    set(gca,'XTick',[-.5 0 .5],'XTickLabel',[]); xlim([-0.5 0.5]); ylabel('PYR','Color',color_pyr);
    
    subplot(5,1,4)
    imagesc_ranked(ts_downstates,[], projectResults.slowOsciResponses.responseZ(sessions_nw,:),[-10 10],...
        projectResults.slowOsciResponses.peakResponseZ(sessions_nw));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    
    plot([win_rebound(1) win_rebound(1)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    plot([win_rebound(2) win_rebound(2)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    set(gca,'XTick',[-.5 0 .5],'XTickLabel',[]); xlim([-0.5 0.5]); ylabel('NW','Color',color_nw);
    
    subplot(5,1,5)
    imagesc_ranked(ts_downstates,[], projectResults.slowOsciResponses.responseZ(sessions_ww,:),[-10 10],...
        projectResults.slowOsciResponses.peakResponseZ(sessions_ww));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    
    plot([win_rebound(1) win_rebound(1)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    plot([win_rebound(2) win_rebound(2)],[0 length(responsive_cells)],'--','color',[.9 .9 .9]);
    xlim([-0.5 0.5]); ylabel('NW','Color',color_nw);
    xlabel('DOWN states (s)');
    colormap jet
    
    % 2.4 DOWN STATES, averages responses figure 6F top
    figure
    plotFill(ts_downstates,projectResults.slowOsciResponses.responseZ(responsive_cells,:),'color',color_cells,'smoothOpt',5);
    plotFill(ts_downstates,projectResults.slowOsciResponses.responseZ(sessions_pyr,:),'color',color_pyr);
    plotFill(ts_downstates,projectResults.slowOsciResponses.responseZ(sessions_nw,:),'color',color_nw);
    plotFill(ts_downstates,projectResults.slowOsciResponses.responseZ(sessions_ww,:),'color',color_ww);
    xlabel('Time since DOWN-state (s)'); ylabel('Rate (s.d.)');
    
    % 
    roundSize = 8;
    figure
    hold on
    rand_pos = ones(length(find(sessions_pyr)),1)+(rand(length(find(sessions_pyr)),1)-0.5)*0.25;
    plot(rand_pos, projectResults.slowOsciResponses.peakResponseZ(sessions_pyr),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    plot(rand_pos + .3, projectResults.slowOsciResponses.peakResponseZRebound(sessions_pyr),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    count = 1;
    for ii = find(sessions_pyr)
        plot([rand_pos(count) rand_pos(count) + .3], [projectResults.slowOsciResponses.peakResponseZ(ii) projectResults.slowOsciResponses.peakResponseZRebound(ii)],'-','LineWidth',.5,'color',[.8 .8 .8]);
        count = 1 + count;
    end
    plot([.9 .9], quantile(projectResults.slowOsciResponses.peakResponseZ(sessions_nw),[.25 .75]),'-','Color',color_pyr);
    plot(.9, median(projectResults.slowOsciResponses.peakResponseZ(sessions_nw)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_pyr,'MarkerFaceColor',color_pyr);
    plot([.9 .9]+.4, quantile(projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw),[.25 .75]),'-','Color',color_pyr);
    plot(.9 + .4, median(projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_pyr,'MarkerFaceColor',[1 1 1]);
    
    rand_pos = ones(length(find(sessions_nw)),1)+(rand(length(find(sessions_nw)),1)-0.5)*0.25;
    plot(rand_pos + 1, projectResults.slowOsciResponses.peakResponseZ(sessions_nw),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    plot(rand_pos + 1.3, projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    count = 1;
    for ii = find(sessions_nw)
        plot([rand_pos(count) rand_pos(count) + .3]+1, [projectResults.slowOsciResponses.peakResponseZ(ii) projectResults.slowOsciResponses.peakResponseZRebound(ii)],'-','LineWidth',.5,'color',[.8 .8 .8]);
        count = 1 + count;
    end
    plot([.9 .9]+1, quantile(projectResults.slowOsciResponses.peakResponseZ(sessions_nw),[.25 .75]),'-','Color',color_nw);
    plot(.9+1, median(projectResults.slowOsciResponses.peakResponseZ(sessions_nw)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_nw,'MarkerFaceColor',color_nw);
    plot([.9 .9] + 1.4, quantile(projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw),[.25 .75]),'-','Color',color_nw);
    plot(.9 + 1.4, median(projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_nw,'MarkerFaceColor',[1 1 1]);
    
    rand_pos = ones(length(find(sessions_ww)),1)+(rand(length(find(sessions_ww)),1)-0.5)*0.25;
    plot(rand_pos + 2, projectResults.slowOsciResponses.peakResponseZ(sessions_ww),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    plot(rand_pos + 2.3, projectResults.slowOsciResponses.peakResponseZRebound(sessions_ww),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    count = 1;
    for ii = find(sessions_ww)
        plot([rand_pos(count) rand_pos(count) + .3]+2, [projectResults.slowOsciResponses.peakResponseZ(ii) projectResults.slowOsciResponses.peakResponseZRebound(ii)],'-','LineWidth',.5,'color',[.8 .8 .8]);
        count = 1 + count;
    end
    plot([.9 .9]+2, quantile(projectResults.slowOsciResponses.peakResponseZ(sessions_ww),[.25 .75]),'-','Color',color_ww);
    plot(.9+2, median(projectResults.slowOsciResponses.peakResponseZ(sessions_ww)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_ww,'MarkerFaceColor',color_ww);
    plot([.9 .9] + 2.4, quantile(projectResults.slowOsciResponses.peakResponseZRebound(sessions_ww),[.25 .75]),'-','Color',color_ww);
    plot(.9 + 2.4, median(projectResults.slowOsciResponses.peakResponseZRebound(sessions_ww)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_ww,'MarkerFaceColor',[1 1 1]);
    
    rand_pos = ones(length(find(responsive_cells)),1)+(rand(length(find(responsive_cells)),1)-0.5)*0.25;
    plot(rand_pos + 3, projectResults.slowOsciResponses.peakResponseZ(responsive_cells),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    plot(rand_pos + 3.3, projectResults.slowOsciResponses.peakResponseZRebound(responsive_cells),'.','MarkerSize',5,'MarkerEdgeColor',[.1 .1 .1],'MarkerFaceColor','none');
    count = 1;
    for ii = find(responsive_cells)
        plot([rand_pos(count) rand_pos(count) + .3]+3, [projectResults.slowOsciResponses.peakResponseZ(ii) projectResults.slowOsciResponses.peakResponseZRebound(ii)],'-','LineWidth',.5,'color',[.8 .8 .8]);
        count = 1 + count;
    end
    plot([.9 .9]+3, quantile(projectResults.slowOsciResponses.peakResponseZ(responsive_cells),[.25 .75]),'-','Color',color_cells);
    plot(.9+3, median(projectResults.slowOsciResponses.peakResponseZ(responsive_cells)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_cells,'MarkerFaceColor',color_cells);
    plot([.9 .9] + 3.4, quantile(projectResults.slowOsciResponses.peakResponseZRebound(responsive_cells),[.25 .75]),'-','Color',color_cells);
    plot(.9 + 3.4, median(projectResults.slowOsciResponses.peakResponseZRebound(responsive_cells)),'o','MarkerSize',roundSize,'MarkerEdgeColor',color_cells,'MarkerFaceColor',[1 1 1]);
    plot([.5 4.7],[0 0],'color',[.7 .7 .7]);
    ylim([-15 12]); xlim([.5 4.7]);
    set(gca,'TickDir','out','XTick',[1:4],'XTickLabel',{'PYR','NW','WW','ID2/DLX'},'XTickLabelRotation',45);
    ylabel('Rate (S.D.)');
    
    groupStats({projectResults.slowOsciResponses.peakResponseZ(sessions_pyr),projectResults.slowOsciResponses.peakResponseZ(sessions_nw),projectResults.slowOsciResponses.peakResponseZ(sessions_ww)...
        ,projectResults.slowOsciResponses.peakResponseZ(responsive_cells)});
    
    groupStats({projectResults.slowOsciResponses.peakResponseZRebound(sessions_pyr),projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw),projectResults.slowOsciResponses.peakResponseZRebound(sessions_ww)...
        ,projectResults.slowOsciResponses.peakResponseZRebound(responsive_cells)});
    
    % 2.5 DOWN-state vs post-DOWN correlation, 6G
    figure
    subplot(4,1,1)
    gc = groupCorr(projectResults.slowOsciResponses.peakResponseZ(sessions_pyr), projectResults.slowOsciResponses.peakResponseZRebound(sessions_pyr),...
        'plotType','XYDispersion','MarkerColor',color_pyr_dark,'inAxis',true,'removeOutliers',true);
    ylabel('post-DOWN (s.d.)');
    set(gca,'XTickLabel',[]);
    xlim([-15 0]);
    
    subplot(4,1,2)
    gc = groupCorr(projectResults.slowOsciResponses.peakResponseZ(sessions_nw), projectResults.slowOsciResponses.peakResponseZRebound(sessions_nw),...
        'plotType','XYDispersion','MarkerColor',color_nw_dark,'inAxis',true,'removeOutliers',true);
    set(gca,'XTickLabel',[]);
    xlim([-15 0]);
    
    subplot(4,1,3)
    gc = groupCorr(projectResults.slowOsciResponses.peakResponseZ(sessions_ww), projectResults.slowOsciResponses.peakResponseZRebound(sessions_ww),...
        'plotType','XYDispersion','MarkerColor',color_ww_dark,'inAxis',true,'removeOutliers',true);
    set(gca,'XTickLabel',[]);
    xlim([-15 0]);
    
    subplot(4,1,4)
    rc = find(responsive_cells); temp = rc(5); rc(5) = rc(9); rc(9) = temp;  temp = rc(6); rc(6) = rc(8); rc(8) = temp;
    gc = groupCorr(projectResults.slowOsciResponses.peakResponseZ(responsive_cells), projectResults.slowOsciResponses.peakResponseZRebound(rc),...
        'plotType','XYDispersion','MarkerColor',color_cells_dark,'inAxis',true,'removeOutliers',false);
    xlim([-15 0]);
    xlabel('DOWN-resp (s.d.)')
    
    
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_nw)), projectResults.optogeneticResponses.latencyAbs(sessions_nw),'plotType','XYDispersion','MarkerColor',color_nw_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10); set(gca,'XTickLabel',[]);
    
    subplot(3,1,3)
    gc = groupCorr(log10(cell_metrics.firingRate(sessions_ww)), projectResults.optogeneticResponses.latencyAbs(sessions_ww),'plotType','XYDispersion','MarkerColor',color_ww_dark,'inAxis',true);
    ylabel('Light response (s.d.)'); xlim(log10([0.1 100]));
    LogScale('x',10);
    

    projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2 = projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime;
    projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2([112   173   175]) = projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime([112   173   175])/2;
    projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2 = projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2 - min(projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2);
    projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2 = projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2/max(projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2);
    
    % 2.6 DOWN to DOWN sequence
    groupStats({projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_pyr),projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_nw),projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_ww)...
        ,projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(responsive_cells)});
    
    figure % coeff of variance
    hold on
    [gs] = groupStats({ projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_pyr),  projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_nw),...
         projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(sessions_ww),  projectResults.slowOsciSpikesRank.median_allEventsMedian_RelTime2(responsive_cells)},[],'color',...
        [color_pyr; color_nw; color_ww; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true,'orientation', 'horizontal','roundPlotSize',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PYR','NW','WW',name_cells},'XTickLabelRotation',45); ylim([.2 .81]); ylabel('UP-DOWN rank order (norm)');

    
    
    
    % 1.4 SLOW GAMMA
    ps = projectResults.lGammaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    for ii = 1:size(projectResults.lGammaModulation.phasedistros,1)
        projectResults.lGammaModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.lGammaModulation.phasedistros(ii,:),[],@mean));
        projectResults.lGammaModulation.phasedistro_lowRes(ii,:) =...
            projectResults.lGammaModulation.phasedistro_lowRes(ii,:)/sum(projectResults.lGammaModulation.phasedistro_lowRes(ii,:));
        projectResults.lGammaModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.lGammaModulation.phasedistro_lowRes(ii,:) projectResults.lGammaModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.lGammaModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.lGammaModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.lGammaModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.lGammaModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.lGammaModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.lGammaModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Theta phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.lGammaModulation.phasestats_m(sessions_pyr), projectResults.lGammaModulation.phasestats_m(sessions_int),...
        projectResults.lGammaModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Theta preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.lGammaModulation.phasestats_r(sessions_pyr)), (projectResults.lGammaModulation.phasestats_r(sessions_int)),...
        (projectResults.lGammaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 FAST GAMMA
    ps = projectResults.lGammaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    for ii = 1:size(projectResults.hGammaModulation.phasedistros,1)
        projectResults.hGammaModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.hGammaModulation.phasedistros(ii,:),[],@mean));
        projectResults.hGammaModulation.phasedistro_lowRes(ii,:) =...
            projectResults.hGammaModulation.phasedistro_lowRes(ii,:)/sum(projectResults.hGammaModulation.phasedistro_lowRes(ii,:));
        projectResults.hGammaModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.hGammaModulation.phasedistro_lowRes(ii,:) projectResults.hGammaModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.hGammaModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.hGammaModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.hGammaModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.hGammaModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.hGammaModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.hGammaModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Theta phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.hGammaModulation.phasestats_m(sessions_pyr), projectResults.hGammaModulation.phasestats_m(sessions_int),...
        projectResults.hGammaModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('GammaM preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.hGammaModulation.phasestats_r(sessions_pyr)), (projectResults.hGammaModulation.phasestats_r(sessions_int)),...
        (projectResults.hGammaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    
%     % 2. DIMENSIONALITY REDUCTION
%     Y = tsne([projectResults.cell_metrics.troughToPeak'...
%         projectResults.cell_metrics.acg_tau_rise'...
%         log10(projectResults.cell_metrics.firingRate_MA)'...
%         projectResults.cell_metrics.firingRateCV'...
%         projectResults.averageCCG.peakResponseZ...
%         projectResults.ripplesResponses.peakResponseZ_norm...
%         projectResults.ripplePhaseModulation.phasestats_m...
%         projectResults.ripplePhaseModulation.phasestats_r...
%         projectResults.thetaModulation.phasestats_m...
%         projectResults.thetaModulation.phasestats_r]);
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),30,color_cells_dark,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     
%     % plot features
%     figure
%     subplot(3,2,1)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7,log10(projectResults.cell_metrics.troughToPeak(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('troughToPeak','FontWeight','normal');
%     subplot(3,2,2)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.acg_tau_rise(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('acg_tau_rise','FontWeight','normal');
%     subplot(3,2,3)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.firingRate_MA(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('firingRate_MA','FontWeight','normal');
%     xlabel('Dimension 2');
%     subplot(3,2,4)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.ripplesResponses.peakResponseZ_norm(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('ripple response','FontWeight','normal');
%     subplot(3,2,5)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.ripplePhaseModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple phase','FontWeight','normal'); xlabel('                                              Dimension 1');
%     subplot(3,2,6)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.thetaModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple response','FontWeight','normal');
    
    
    
    
%     % which are the other cells?
%     responsive_cells_allSessions = any(projectResults.optogeneticResponses.threeWaysTest'==1);
%     
%     targetSessCells_ivy = strcmpi(projectResults.geneticLine,'id2/nkx2.1/ai80');
%     targetSessCells_pv = strcmpi(projectResults.geneticLine,'pv/ai32');
%     targetSessCells_sst = strcmpi(projectResults.geneticLine,'sst/ai32');
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),50,color_cells_dark,'LineWidth',2);
%     
%     scatter(Y(responsive_cells_allSessions(targetSessCells_ivy),1),Y(responsive_cells_allSessions(targetSessCells_ivy),2),50,color_id2nkx,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     scatter(Y(responsive_cells_allSessions(targetSessCells_pv),1),Y(responsive_cells_allSessions(targetSessCells_pv),2),50,color_pv,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
% %     scatter(Y(responsive_cells_allSessions(targetSessCells_sst),1),Y(responsive_cells_allSessions(targetSessCells_sst),2),50,color_sst,'LineWidth',2);
% %     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
end

% CHAPTER_2: PARVALBUMIN
for z = 1
    % PV cells features
    targetSessCells = strcmpi(projectResults.geneticLine,'pv/ai32') & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    allTargetSessCells.pv = targetSessCells;
    
    name_cells = 'PV+';
    color_cells = color_pv;
    color_cells_dark = color_pv_dark;
    
    % 1.1 Light responses
    ts = projectResults.optogeneticResponses.timestamps;
    responsive_cells = any(projectResults.optogeneticResponses.threeWaysTest'==1) & targetSessCells & (projectResults.optogeneticResponses.checkedCells==1)';
    responsive_cells([1903 2039 2579]) = 0;
    
    allResponsive_cells.pv = responsive_cells;
    [~,optimal_pulse] = max(projectResults.optogeneticResponses.rateDuringPulse');
    sessions_pyr = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_pyr;
    sessions_int = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_int;
    all_pyr.pv = sessions_pyr;
    all_int.pv = sessions_int;
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(responsive_cells,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt));
    xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
    xlabel('Time since light stimulation (ms)'); colormap jet
    
    % 
    [gs] = groupStats({projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt), projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt),...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Light responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.2 INTRINSIC FEATURES
    % A) SPIKE WIDTH
    cell_metrics = projectResults.cell_metrics;
    waveforms_timestmaps = cell_metrics.waveforms.time{1};
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');

    subplot(1,2,2)
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(sessions_pyr), projectResults.cell_metrics.troughToPeak(sessions_int),...
        projectResults.cell_metrics.troughToPeak(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % B) ACG
    acg_timestmaps = linspace(-50,50, size(cell_metrics.acg.narrow_normalized,1));
    for ii = 1:size(cell_metrics.acg.narrow_normalized,2)
        cell_metrics.acg.narrow_probability(:,ii) = smooth(cell_metrics.acg.narrow(:,ii)/sum(cell_metrics.acg.narrow(:,ii)),10);
    end
    cell_metrics.acg.narrow_probability([100:102],:) = 0;
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    
    subplot(1,2,2)
    [gs] = groupStats({cell_metrics.acg_tau_rise(sessions_pyr), cell_metrics.acg_tau_rise(sessions_int),...
        cell_metrics.acg_tau_rise(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Tau rise (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % log10ACG
    acg_timestmaps = projectResults.acgPeak.acg_time(1,:);
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_pyr,:),1),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_int,:),1),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(responsive_cells,:),1),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    xlim(log10([0.001 1])); LogScale('x',10); xlim(log10([0.0015 1]));
    
    subplot(1,2,2)
    [gs] = groupStats({projectResults.acgPeak.peakTime(sessions_pyr), projectResults.acgPeak.peakTime(sessions_int),...
        projectResults.acgPeak.peakTime(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylim(log10([0.001 1])); LogScale('y',10); ylim(log10([0.0015 1]));
    ylabel('ACG peak (s)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % C) FIRING RATE PER STATES
    figure
    subplot(1,2,1)
    [gs] = groupStats({log10(cell_metrics.firingRate_MA(sessions_pyr)), log10(cell_metrics.firingRate_MA(sessions_int)),...
        log10(cell_metrics.firingRate_MA(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Firing rate (Hz)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.01 100])); LogScale('y',10);
    
    subplot(1,2,2)
    [gs] = groupStats({log10(cell_metrics.firingRateCV(sessions_pyr)), log10(cell_metrics.firingRateCV(sessions_int)),...
        log10(cell_metrics.firingRateCV(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Coefficient of variance');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.1 10])); LogScale('y',10);
    
    figure
    subplot(1,3,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot([1 2 3 4]-0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_pyr); p.Color(4) = .05;
    end
    
    for ii = find(sessions_int)
        p = plot([1 2 3 4], [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_int); p.Color(4) = .05;
    end
    
    for ii = find(responsive_cells)
        p = plot([1 2 3 4]+0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_cells); p.Color(4) = .05;
    end
    
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_pyr)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_pyr)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_pyr)), log10(cell_metrics.firingRate_REMstate(sessions_pyr))},[],'color',[color_pyr],...
        'plotData',true,'plotType','medianBall','labelSummary',false,'inAxis',true,'posOffset',[-0.1],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_int)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_int)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_int)), log10(cell_metrics.firingRate_REMstate(sessions_int))},[],'color',[color_int],...
        'plotData',true,'plotType','medianBall','labelSummary',false,'inAxis',true,'posOffset',[0],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(responsive_cells)), log10(cell_metrics.firingRate_WAKEnontheta(responsive_cells)),...
        log10(cell_metrics.firingRate_NREMstate(responsive_cells)), log10(cell_metrics.firingRate_REMstate(responsive_cells))},[],'color',[color_cells],...
        'plotData',true,'plotType','medianBall','labelSummary',false,'inAxis',true,'posOffset',[0.1],'sigStar',false);
    ylim(log10([0.05 40])); LogScale('y',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Run','Quiet','NREM','REM'},'XTickLabelRotation',45);
    ylim(log10([0.1 15])); LogScale('y',10);
    ylabel('Firing rate (Hz)');
    
    cell_metrics.run_quiet = (cell_metrics.firingRate_WAKEtheta - cell_metrics.firingRate_WAKEnontheta)./(cell_metrics.firingRate_WAKEtheta + cell_metrics.firingRate_WAKEnontheta);
    cell_metrics.run_quiet(cell_metrics.run_quiet==-1 | cell_metrics.run_quiet==1) = NaN;
    
    cell_metrics.rem_nrem = (-cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate)./(cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate);
    cell_metrics.rem_nrem(cell_metrics.rem_nrem==-1 | cell_metrics.rem_nrem==1) = NaN;
    
    subplot(1,3,2)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.run_quiet(sessions_pyr), cell_metrics.run_quiet(sessions_int),...
        cell_metrics.run_quiet(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Quiet <----> Run'); ylim([-1 1]);
    
    subplot(1,3,3)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.rem_nrem(sessions_pyr), cell_metrics.rem_nrem(sessions_int),...
        cell_metrics.rem_nrem(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('NREM <----> REM'); ylim([-1 1]);
    
    % 1.3 AVERAGE CCG
    win_resp = [-0.03 0.03];
    ts_CCG = projectResults.averageCCG.timestamps;
    
    win_Z = find(ts_CCG<=-0.1);
    for ii = 1:size(projectResults.averageCCG.ZmeanCCG,1)
        projectResults.averageCCG.ZmeanCCGSmooth(ii,:) = smooth(projectResults.averageCCG.ZmeanCCG(ii,:),3);
    end
    
    win = find(ts_CCG>=win_resp(1) & ts_CCG<=win_resp(2));
    projectResults.averageCCG.peakResponse = nanmean(projectResults.averageCCG.meanCCG(:,win),2); % delta peak response
    projectResults.averageCCG.peakResponseZ = nanmean(projectResults.averageCCG.ZmeanCCG(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(responsive_cells,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_pyr,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_int,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('CCG population responses (ms)');
    colormap jet
    
    [gs] = groupStats({projectResults.averageCCG.peakResponseZ(sessions_pyr), projectResults.averageCCG.peakResponseZ(sessions_int),...
        projectResults.averageCCG.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Average CCG resp (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    % 1.4 RIPPLES FIRING
    win_resp = [-0.025 0.025];
    ts_ripples = projectResults.ripplesResponses.timestamps;
    
    win_Z = find(ts_ripples<=-0.1);
    for ii = 1:size(projectResults.ripplesResponses.responsecurveSmooth,1)
        projectResults.ripplesResponses.responseZ(ii,:) = (projectResults.ripplesResponses.responsecurveSmooth(ii,:) - ...
            mean(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z)))./std(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z));
    end
    
    win = find(ts_ripples>=win_resp(1) & ts_ripples<=win_resp(2));
    projectResults.ripplesResponses.peakResponse = nanmean(projectResults.ripplesResponses.responsecurve(:,win),2); % delta peak response
    projectResults.ripplesResponses.peakResponseZ = nanmean(projectResults.ripplesResponses.responseZ(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(responsive_cells,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_pyr,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_int,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('Ripple responses (ms)');
    colormap jet
    
    % 
    [gs] = groupStats({projectResults.ripplesResponses.peakResponseZ(sessions_pyr), projectResults.ripplesResponses.peakResponseZ(sessions_int),...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Ripple responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 RIPPLES PHASE
    ps = projectResults.ripplePhaseModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    projectResults.ripplePhaseModulation.phases_lowRes = ps_lowRes;
    
    for ii = 1:size(projectResults.ripplePhaseModulation.phasedistros,1)
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.ripplePhaseModulation.phasedistros(ii,:),[],@mean));
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)/sum(projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:));
        projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Ripple phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr), projectResults.ripplePhaseModulation.phasestats_m(sessions_int),...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true,'sigStar',false);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Ripple phase preference (rad)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.ripplePhaseModulation.phasestats_r(sessions_pyr)), (projectResults.ripplePhaseModulation.phasestats_r(sessions_int)),...
        (projectResults.ripplePhaseModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 THETA
    ps = projectResults.thetaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    for ii = 1:size(projectResults.thetaModulation.phasedistros,1)
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.thetaModulation.phasedistros(ii,:),[],@mean));
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            projectResults.thetaModulation.phasedistro_lowRes(ii,:)/sum(projectResults.thetaModulation.phasedistro_lowRes(ii,:));
        projectResults.thetaModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.thetaModulation.phasedistro_lowRes(ii,:) projectResults.thetaModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Theta phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.thetaModulation.phasestats_m(sessions_pyr), projectResults.thetaModulation.phasestats_m(sessions_int),...
        projectResults.thetaModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Theta preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.thetaModulation.phasestats_r(sessions_pyr)), (projectResults.thetaModulation.phasestats_r(sessions_int)),...
        (projectResults.thetaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    
%     % 2. DIMENSIONALITY REDUCTION
%     Y = tsne([projectResults.cell_metrics.troughToPeak'...
%         projectResults.cell_metrics.acg_tau_rise'...
%         log10(projectResults.cell_metrics.firingRate_MA)'...
%         projectResults.cell_metrics.firingRateCV'...
%         projectResults.averageCCG.peakResponseZ...
%         projectResults.ripplesResponses.peakResponseZ_norm...
%         projectResults.ripplePhaseModulation.phasestats_m...
%         projectResults.ripplePhaseModulation.phasestats_r...
%         projectResults.thetaModulation.phasestats_m...
%         projectResults.thetaModulation.phasestats_r]);
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),30,color_cells_dark,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     
%     % plot features
%     figure
%     subplot(3,2,1)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7,log10(projectResults.cell_metrics.troughToPeak(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('troughToPeak','FontWeight','normal');
%     subplot(3,2,2)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.acg_tau_rise(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('acg_tau_rise','FontWeight','normal');
%     subplot(3,2,3)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.firingRate_MA(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('firingRate_MA','FontWeight','normal');
%     xlabel('Dimension 2');
%     subplot(3,2,4)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.ripplesResponses.peakResponseZ_norm(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('ripple response','FontWeight','normal');
%     subplot(3,2,5)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.ripplePhaseModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple phase','FontWeight','normal'); xlabel('                                              Dimension 1');
%     subplot(3,2,6)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.thetaModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple response','FontWeight','normal');
    
    
    
    
%     % which are the other cells?
%     responsive_cells_allSessions = any(projectResults.optogeneticResponses.threeWaysTest'==1);
%     
%     targetSessCells_ivy = strcmpi(projectResults.geneticLine,'id2/nkx2.1/ai80');
%     targetSessCells_pv = strcmpi(projectResults.geneticLine,'pv/ai32');
%     targetSessCells_sst = strcmpi(projectResults.geneticLine,'sst/ai32');
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),50,color_cells_dark,'LineWidth',2);
%     
%     scatter(Y(responsive_cells_allSessions(targetSessCells_ivy),1),Y(responsive_cells_allSessions(targetSessCells_ivy),2),50,color_id2nkx,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     scatter(Y(responsive_cells_allSessions(targetSessCells_pv),1),Y(responsive_cells_allSessions(targetSessCells_pv),2),50,color_pv,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
% %     scatter(Y(responsive_cells_allSessions(targetSessCells_sst),1),Y(responsive_cells_allSessions(targetSessCells_sst),2),50,color_sst,'LineWidth',2);
% %     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
end

% CHAPTER_3: SST
for z = 1
    % SST cells features
    targetSessCells = strcmpi(projectResults.geneticLine,'sst/ai32') & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    allTargetSessCells.sst = targetSessCells;
    
    name_cells = 'SST+';
    color_cells = color_sst;
    color_cells_dark = color_sst_dark;
    
    % 1.1 Light responses
    ts = projectResults.optogeneticResponses.timestamps;
    responsive_cells = any(projectResults.optogeneticResponses.threeWaysTest'==1) & targetSessCells & (projectResults.optogeneticResponses.checkedCells==1)';
    % responsive_cells([1903 2039 2579]) = 0;
    
    allResponsive_cells.sst = responsive_cells;
    [~,optimal_pulse] = max(projectResults.optogeneticResponses.rateDuringPulse');
    sessions_pyr = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_pyr;
    sessions_int = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_int;
    all_pyr.sst = sessions_pyr;
    all_int.sst = sessions_int;
    
    figure % 2.5 x 7
    post_opt_local = 3;
    subplot(4,1,[1 2])
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(responsive_cells,post_opt_local,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt_local));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt_local,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt_local));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt_local,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt_local));
    xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
    xlabel('Time since light stimulation (ms)'); colormap jet
    
    % 
    [gs] = groupStats({projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt), projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt),...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Light responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.2 INTRINSIC FEATURES
    % A) SPIKE WIDTH
    cell_metrics = projectResults.cell_metrics;
    waveforms_timestmaps = cell_metrics.waveforms.time{1};
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');

    subplot(1,2,2)
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(sessions_pyr), projectResults.cell_metrics.troughToPeak(sessions_int),...
        projectResults.cell_metrics.troughToPeak(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % B) ACG
    acg_timestmaps = linspace(-50,50, size(cell_metrics.acg.narrow_normalized,1));
    for ii = 1:size(cell_metrics.acg.narrow_normalized,2)
        cell_metrics.acg.narrow_probability(:,ii) = smooth(cell_metrics.acg.narrow(:,ii)/sum(cell_metrics.acg.narrow(:,ii)),10);
    end
    cell_metrics.acg.narrow_probability([100:102],:) = 0;
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    
    subplot(1,2,2)
    [gs] = groupStats({cell_metrics.acg_tau_rise(sessions_pyr), cell_metrics.acg_tau_rise(sessions_int),...
        cell_metrics.acg_tau_rise(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Tau rise (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % log10ACG
    acg_timestmaps = projectResults.acgPeak.acg_time(1,:);
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_pyr,:),1),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_int,:),1),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(responsive_cells,:),1),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    xlim(log10([0.001 1])); LogScale('x',10); xlim(log10([0.0015 1]));
    
    subplot(1,2,2)
    [gs] = groupStats({projectResults.acgPeak.peakTime(sessions_pyr), projectResults.acgPeak.peakTime(sessions_int),...
        projectResults.acgPeak.peakTime(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylim(log10([0.001 1])); LogScale('y',10); ylim(log10([0.0015 1]));
    ylabel('ACG peak (s)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % C) FIRING RATE PER STATES
    figure
    subplot(1,2,1)
    [gs] = groupStats({log10(cell_metrics.firingRate_MA(sessions_pyr)), log10(cell_metrics.firingRate_MA(sessions_int)),...
        log10(cell_metrics.firingRate_MA(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Firing rate (Hz)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.01 100])); LogScale('y',10);
    
    subplot(1,2,2)
    [gs] = groupStats({log10(cell_metrics.firingRateCV(sessions_pyr)), log10(cell_metrics.firingRateCV(sessions_int)),...
        log10(cell_metrics.firingRateCV(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Coefficient of variance');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.1 10])); LogScale('y',10);
    
    figure
    subplot(1,3,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot([1 2 3 4]-0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_pyr); p.Color(4) = .05;
    end
    
    for ii = find(sessions_int)
        p = plot([1 2 3 4], [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_int); p.Color(4) = .05;
    end
    
    for ii = find(responsive_cells)
        p = plot([1 2 3 4]+0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_cells); p.Color(4) = .05;
    end
    
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_pyr)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_pyr)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_pyr)), log10(cell_metrics.firingRate_REMstate(sessions_pyr))},[],'color',[color_pyr],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[-0.1],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_int)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_int)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_int)), log10(cell_metrics.firingRate_REMstate(sessions_int))},[],'color',[color_int],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(responsive_cells)), log10(cell_metrics.firingRate_WAKEnontheta(responsive_cells)),...
        log10(cell_metrics.firingRate_NREMstate(responsive_cells)), log10(cell_metrics.firingRate_REMstate(responsive_cells))},[],'color',[color_cells],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0.1],'sigStar',false);
    ylim(log10([0.05 40])); LogScale('y',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Run','Quiet','NREM','REM'},'XTickLabelRotation',45);
    ylim(log10([0.1 15])); LogScale('y',10);
    ylabel('Firing rate (Hz)');
    
    cell_metrics.run_quiet = (cell_metrics.firingRate_WAKEtheta - cell_metrics.firingRate_WAKEnontheta)./(cell_metrics.firingRate_WAKEtheta + cell_metrics.firingRate_WAKEnontheta);
    cell_metrics.run_quiet(cell_metrics.run_quiet==-1 | cell_metrics.run_quiet==1) = NaN;
    
    cell_metrics.rem_nrem = (-cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate)./(cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate);
    cell_metrics.rem_nrem(cell_metrics.rem_nrem==-1 | cell_metrics.rem_nrem==1) = NaN;
    
    subplot(1,3,2)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.run_quiet(sessions_pyr), cell_metrics.run_quiet(sessions_int),...
        cell_metrics.run_quiet(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Quiet <----> Run'); ylim([-1 1]);
    
    subplot(1,3,3)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.rem_nrem(sessions_pyr), cell_metrics.rem_nrem(sessions_int),...
        cell_metrics.rem_nrem(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('NREM <----> REM'); ylim([-1 1]);
    
    % 1.3 AVERAGE CCG
    win_resp = [-0.03 0.03];
    ts_CCG = projectResults.averageCCG.timestamps;
    
    win_Z = find(ts_CCG<=-0.1);
    for ii = 1:size(projectResults.averageCCG.ZmeanCCG,1)
        projectResults.averageCCG.ZmeanCCGSmooth(ii,:) = smooth(projectResults.averageCCG.ZmeanCCG(ii,:),3);
    end
    
    win = find(ts_CCG>=win_resp(1) & ts_CCG<=win_resp(2));
    projectResults.averageCCG.peakResponse = nanmean(projectResults.averageCCG.meanCCG(:,win),2); % delta peak response
    projectResults.averageCCG.peakResponseZ = nanmean(projectResults.averageCCG.ZmeanCCG(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(responsive_cells,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_pyr,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_int,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('CCG population responses (ms)');
    colormap jet
    
    [gs] = groupStats({projectResults.averageCCG.peakResponseZ(sessions_pyr), projectResults.averageCCG.peakResponseZ(sessions_int),...
        projectResults.averageCCG.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Average CCG resp (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    % 1.4 RIPPLES FIRING
    win_resp = [-0.025 0.025];
    ts_ripples = projectResults.ripplesResponses.timestamps;
    
    win_Z = find(ts_ripples<=-0.1);
    for ii = 1:size(projectResults.ripplesResponses.responsecurveSmooth,1)
        projectResults.ripplesResponses.responseZ(ii,:) = (projectResults.ripplesResponses.responsecurveSmooth(ii,:) - ...
            mean(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z)))./std(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z));
    end
    
    win = find(ts_ripples>=win_resp(1) & ts_ripples<=win_resp(2));
    projectResults.ripplesResponses.peakResponse = nanmean(projectResults.ripplesResponses.responsecurve(:,win),2); % delta peak response
    projectResults.ripplesResponses.peakResponseZ = nanmean(projectResults.ripplesResponses.responseZ(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(responsive_cells,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_pyr,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_int,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('Ripple responses (ms)');
    colormap jet
    
    % 
    [gs] = groupStats({projectResults.ripplesResponses.peakResponseZ(sessions_pyr), projectResults.ripplesResponses.peakResponseZ(sessions_int),...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Ripple responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 RIPPLES PHASE
    ps = projectResults.ripplePhaseModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    projectResults.ripplePhaseModulation.phases_lowRes = ps_lowRes;
    
    for ii = 1:size(projectResults.ripplePhaseModulation.phasedistros,1)
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.ripplePhaseModulation.phasedistros(ii,:),[],@mean));
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)/sum(projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:));
        projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Ripple phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr), projectResults.ripplePhaseModulation.phasestats_m(sessions_int),...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true,'sigStar',false);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Ripple phase preference (rad)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.ripplePhaseModulation.phasestats_r(sessions_pyr)), (projectResults.ripplePhaseModulation.phasestats_r(sessions_int)),...
        (projectResults.ripplePhaseModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 THETA
    ps = projectResults.thetaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    for ii = 1:size(projectResults.thetaModulation.phasedistros,1)
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.thetaModulation.phasedistros(ii,:),[],@mean));
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            projectResults.thetaModulation.phasedistro_lowRes(ii,:)/sum(projectResults.thetaModulation.phasedistro_lowRes(ii,:));
        projectResults.thetaModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.thetaModulation.phasedistro_lowRes(ii,:) projectResults.thetaModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Theta phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.thetaModulation.phasestats_m(sessions_pyr), projectResults.thetaModulation.phasestats_m(sessions_int),...
        projectResults.thetaModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Theta preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.thetaModulation.phasestats_r(sessions_pyr)), (projectResults.thetaModulation.phasestats_r(sessions_int)),...
        (projectResults.thetaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    
%     % 2. DIMENSIONALITY REDUCTION
%     Y = tsne([projectResults.cell_metrics.troughToPeak'...
%         projectResults.cell_metrics.acg_tau_rise'...
%         log10(projectResults.cell_metrics.firingRate_MA)'...
%         projectResults.cell_metrics.firingRateCV'...
%         projectResults.averageCCG.peakResponseZ...
%         projectResults.ripplesResponses.peakResponseZ_norm...
%         projectResults.ripplePhaseModulation.phasestats_m...
%         projectResults.ripplePhaseModulation.phasestats_r...
%         projectResults.thetaModulation.phasestats_m...
%         projectResults.thetaModulation.phasestats_r]);
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),30,color_cells_dark,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     
%     % plot features
%     figure
%     subplot(3,2,1)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7,log10(projectResults.cell_metrics.troughToPeak(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('troughToPeak','FontWeight','normal');
%     subplot(3,2,2)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.acg_tau_rise(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('acg_tau_rise','FontWeight','normal');
%     subplot(3,2,3)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.firingRate_MA(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('firingRate_MA','FontWeight','normal');
%     xlabel('Dimension 2');
%     subplot(3,2,4)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.ripplesResponses.peakResponseZ_norm(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('ripple response','FontWeight','normal');
%     subplot(3,2,5)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.ripplePhaseModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple phase','FontWeight','normal'); xlabel('                                              Dimension 1');
%     subplot(3,2,6)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.thetaModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple response','FontWeight','normal');
    
    
    
    
%     % which are the other cells?
%     responsive_cells_allSessions = any(projectResults.optogeneticResponses.threeWaysTest'==1);
%     
%     targetSessCells_ivy = strcmpi(projectResults.geneticLine,'id2/nkx2.1/ai80');
%     targetSessCells_pv = strcmpi(projectResults.geneticLine,'pv/ai32');
%     targetSessCells_sst = strcmpi(projectResults.geneticLine,'sst/ai32');
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),50,color_cells_dark,'LineWidth',2);
%     
%     scatter(Y(responsive_cells_allSessions(targetSessCells_ivy),1),Y(responsive_cells_allSessions(targetSessCells_ivy),2),50,color_id2nkx,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     scatter(Y(responsive_cells_allSessions(targetSessCells_pv),1),Y(responsive_cells_allSessions(targetSessCells_pv),2),50,color_pv,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
% %     scatter(Y(responsive_cells_allSessions(targetSessCells_sst),1),Y(responsive_cells_allSessions(targetSessCells_sst),2),50,color_sst,'LineWidth',2);
% %     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
end

% CHAPTER_4: CAMK2
for z = 1
    % PYR cells features
    targetSessCells = strcmpi(projectResults.geneticLine,'camkii/ai32') & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    allTargetSessCells.camk2 = targetSessCells;
    
    name_cells = 'CaMKII+';
    color_cells = color_camk2;
    color_cells_dark = color_camk2_dark;
    
    % 1.1 Light responses
    ts = projectResults.optogeneticResponses.timestamps;
    responsive_cells = any(projectResults.optogeneticResponses.threeWaysTest'==1) & targetSessCells & (projectResults.optogeneticResponses.checkedCells==1)';
    responsive_cells([2635        2639        2642        2645]) = 0;
    
    allResponsive_cells.camk2 = responsive_cells;
    [~,optimal_pulse] = max(projectResults.optogeneticResponses.rateDuringPulse');
    sessions_pyr = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_pyr;
    sessions_int = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_int;
    all_pyr.camk2 = sessions_pyr;
    all_int.camk2 = sessions_int;
    
    figure % 2.5 x 7
    post_opt_local = 3;
    subplot(4,1,[1 2])
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(responsive_cells,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_pyr,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt));
    set(gca,'XTick',[]); xlim([-0.1 0.5]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts,[], squeeze(nanmean(projectResults.optogeneticResponses.responsecurveZSmooth(sessions_int,post_opt,:),2)),[-5 5],...
        projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt));
    xlim([-0.1 0.5]); ylabel('INT','Color',color_int);
    xlabel('Time since light stimulation (ms)'); colormap jet
    
    % 
    [gs] = groupStats({projectResults.optogeneticResponses.rateZDuringPulse(sessions_pyr,post_opt), projectResults.optogeneticResponses.rateZDuringPulse(sessions_int,post_opt),...
        projectResults.optogeneticResponses.rateZDuringPulse(responsive_cells,post_opt)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Light responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.2 INTRINSIC FEATURES
    % A) SPIKE WIDTH
    cell_metrics = projectResults.cell_metrics;
    waveforms_timestmaps = cell_metrics.waveforms.time{1};
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(waveforms_timestmaps, mean(cell_metrics.waveforms.filt_zscored(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');

    subplot(1,2,2)
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(sessions_pyr), projectResults.cell_metrics.troughToPeak(sessions_int),...
        projectResults.cell_metrics.troughToPeak(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % B) ACG
    acg_timestmaps = linspace(-50,50, size(cell_metrics.acg.narrow_normalized,1));
    for ii = 1:size(cell_metrics.acg.narrow_normalized,2)
        cell_metrics.acg.narrow_probability(:,ii) = smooth(cell_metrics.acg.narrow(:,ii)/sum(cell_metrics.acg.narrow(:,ii)),10);
    end
    cell_metrics.acg.narrow_probability([100:102],:) = 0;
    
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, cell_metrics.acg.narrow_probability(:,ii),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_pyr),2),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,sessions_int),2),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(cell_metrics.acg.narrow_probability(:,responsive_cells),2),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    
    subplot(1,2,2)
    [gs] = groupStats({cell_metrics.acg_tau_rise(sessions_pyr), cell_metrics.acg_tau_rise(sessions_int),...
        cell_metrics.acg_tau_rise(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Tau rise (ms)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % log10ACG
    acg_timestmaps = projectResults.acgPeak.acg_time(1,:);
    figure
    subplot(1,2,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_pyr); p.Color(4) = .05;
    end
    for ii = find(sessions_int)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_int); p.Color(4) = .05;
    end
    for ii = find(responsive_cells)
        p = plot(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(ii,:),'color',color_cells); p.Color(4) = 1;
    end
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_pyr,:),1),'color',color_pyr_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(sessions_int,:),1),'color',color_int_dark,'LineWidth',2);
    plot(acg_timestmaps, mean(projectResults.acgPeak.acg_smoothed_norm(responsive_cells,:),1),'color',color_cells_dark,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Prob'); set(gca, 'TickDir', 'out');
    xlim(log10([0.001 1])); LogScale('x',10); xlim(log10([0.0015 1]));
    
    subplot(1,2,2)
    [gs] = groupStats({projectResults.acgPeak.peakTime(sessions_pyr), projectResults.acgPeak.peakTime(sessions_int),...
        projectResults.acgPeak.peakTime(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylim(log10([0.001 1])); LogScale('y',10); ylim(log10([0.0015 1]));
    ylabel('ACG peak (s)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % C) FIRING RATE PER STATES
    figure
    subplot(1,2,1)
    [gs] = groupStats({log10(cell_metrics.firingRate_MA(sessions_pyr)), log10(cell_metrics.firingRate_MA(sessions_int)),...
        log10(cell_metrics.firingRate_MA(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Firing rate (Hz)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.01 100])); LogScale('y',10);
    
    subplot(1,2,2)
    [gs] = groupStats({log10(cell_metrics.firingRateCV(sessions_pyr)), log10(cell_metrics.firingRateCV(sessions_int)),...
        log10(cell_metrics.firingRateCV(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Coefficient of variance');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    ylim(log10([0.1 10])); LogScale('y',10);
    
    figure
    subplot(1,3,1)
    hold on
    for ii = find(sessions_pyr)
        p = plot([1 2 3 4]-0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_pyr); p.Color(4) = .05;
    end
    
    for ii = find(sessions_int)
        p = plot([1 2 3 4], [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_int); p.Color(4) = .05;
    end
    
    for ii = find(responsive_cells)
        p = plot([1 2 3 4]+0.1, [log10(cell_metrics.firingRate_WAKEtheta(ii)), log10(cell_metrics.firingRate_WAKEnontheta(ii)),...
            log10(cell_metrics.firingRate_NREMstate(ii)),log10(cell_metrics.firingRate_REMstate(ii))],'color',color_cells); p.Color(4) = .05;
    end
    
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_pyr)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_pyr)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_pyr)), log10(cell_metrics.firingRate_REMstate(sessions_pyr))},[],'color',[color_pyr],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[-0.1],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_int)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_int)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_int)), log10(cell_metrics.firingRate_REMstate(sessions_int))},[],'color',[color_int],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(responsive_cells)), log10(cell_metrics.firingRate_WAKEnontheta(responsive_cells)),...
        log10(cell_metrics.firingRate_NREMstate(responsive_cells)), log10(cell_metrics.firingRate_REMstate(responsive_cells))},[],'color',[color_cells],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0.1],'sigStar',false);
    ylim(log10([0.05 40])); LogScale('y',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Run','Quiet','NREM','REM'},'XTickLabelRotation',45);
    ylim(log10([0.1 15])); LogScale('y',10);
    ylabel('Firing rate (Hz)');
    
    cell_metrics.run_quiet = (cell_metrics.firingRate_WAKEtheta - cell_metrics.firingRate_WAKEnontheta)./(cell_metrics.firingRate_WAKEtheta + cell_metrics.firingRate_WAKEnontheta);
    cell_metrics.run_quiet(cell_metrics.run_quiet==-1 | cell_metrics.run_quiet==1) = NaN;
    
    cell_metrics.rem_nrem = (-cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate)./(cell_metrics.firingRate_NREMepisode + cell_metrics.firingRate_REMstate);
    cell_metrics.rem_nrem(cell_metrics.rem_nrem==-1 | cell_metrics.rem_nrem==1) = NaN;
    
    subplot(1,3,2)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.run_quiet(sessions_pyr), cell_metrics.run_quiet(sessions_int),...
        cell_metrics.run_quiet(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Quiet <----> Run'); ylim([-1 1]);
    
    subplot(1,3,3)
    hold on
    plot([0 4],[0 0],'color',[.7 .7 .7]);
    [gs] = groupStats({cell_metrics.rem_nrem(sessions_pyr), cell_metrics.rem_nrem(sessions_int),...
        cell_metrics.rem_nrem(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('NREM <----> REM'); ylim([-1 1]);
    
    % 1.3 AVERAGE CCG
    win_resp = [-0.03 0.03];
    ts_CCG = projectResults.averageCCG.timestamps;
    
    win_Z = find(ts_CCG<=-0.1);
    for ii = 1:size(projectResults.averageCCG.ZmeanCCG,1)
        projectResults.averageCCG.ZmeanCCGSmooth(ii,:) = smooth(projectResults.averageCCG.ZmeanCCG(ii,:),3);
    end
    
    win = find(ts_CCG>=win_resp(1) & ts_CCG<=win_resp(2));
    projectResults.averageCCG.peakResponse = nanmean(projectResults.averageCCG.meanCCG(:,win),2); % delta peak response
    projectResults.averageCCG.peakResponseZ = nanmean(projectResults.averageCCG.ZmeanCCG(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(responsive_cells,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_pyr,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_CCG,[], projectResults.averageCCG.ZmeanCCGSmooth(sessions_int,:),[-3 3],...
        projectResults.averageCCG.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('CCG population responses (ms)');
    colormap jet
    
    [gs] = groupStats({projectResults.averageCCG.peakResponseZ(sessions_pyr), projectResults.averageCCG.peakResponseZ(sessions_int),...
        projectResults.averageCCG.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Average CCG resp (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    % 1.4 RIPPLES FIRING
    win_resp = [-0.025 0.025];
    ts_ripples = projectResults.ripplesResponses.timestamps;
    
    win_Z = find(ts_ripples<=-0.1);
    for ii = 1:size(projectResults.ripplesResponses.responsecurveSmooth,1)
        projectResults.ripplesResponses.responseZ(ii,:) = (projectResults.ripplesResponses.responsecurveSmooth(ii,:) - ...
            mean(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z)))./std(projectResults.ripplesResponses.responsecurveSmooth(ii,win_Z));
    end
    
    win = find(ts_ripples>=win_resp(1) & ts_ripples<=win_resp(2));
    projectResults.ripplesResponses.peakResponse = nanmean(projectResults.ripplesResponses.responsecurve(:,win),2); % delta peak response
    projectResults.ripplesResponses.peakResponseZ = nanmean(projectResults.ripplesResponses.responseZ(:,win),2); % delta peak response
    
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(responsive_cells,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_pyr,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_pyr));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    set(gca,'XTick',[]); xlim([-0.3 0.3]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked(ts_ripples,[], projectResults.ripplesResponses.responseZ(sessions_int,:),[-20 20],...
        projectResults.ripplesResponses.peakResponseZ(sessions_int));
    hold on
    plot([win_resp(1) win_resp(1)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 length(responsive_cells)],'color',[.9 .9 .9]);
    xlim([-0.3 0.3]); ylabel('INT','Color',color_int);
    xlabel('Ripple responses (ms)');
    colormap jet
    
    % 
    [gs] = groupStats({projectResults.ripplesResponses.peakResponseZ(sessions_pyr), projectResults.ripplesResponses.peakResponseZ(sessions_int),...
        projectResults.ripplesResponses.peakResponseZ(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Ripple responses (SD)');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 RIPPLES PHASE
    ps = projectResults.ripplePhaseModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    projectResults.ripplePhaseModulation.phases_lowRes = ps_lowRes;
    
    for ii = 1:size(projectResults.ripplePhaseModulation.phasedistros,1)
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.ripplePhaseModulation.phasedistros(ii,:),[],@mean));
        projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) =...
            projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)/sum(projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:));
        projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:) projectResults.ripplePhaseModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.ripplePhaseModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Ripple phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.ripplePhaseModulation.phasestats_m(sessions_pyr), projectResults.ripplePhaseModulation.phasestats_m(sessions_int),...
        projectResults.ripplePhaseModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true,'sigStar',false);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Ripple phase preference (rad)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.ripplePhaseModulation.phasestats_r(sessions_pyr)), (projectResults.ripplePhaseModulation.phasestats_r(sessions_int)),...
        (projectResults.ripplePhaseModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    % 1.4 THETA
    ps = projectResults.thetaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    for ii = 1:size(projectResults.thetaModulation.phasedistros,1)
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            wrapTo2Pi(accumarray(ps_bins',projectResults.thetaModulation.phasedistros(ii,:),[],@mean));
        projectResults.thetaModulation.phasedistro_lowRes(ii,:) =...
            projectResults.thetaModulation.phasedistro_lowRes(ii,:)/sum(projectResults.thetaModulation.phasedistro_lowRes(ii,:));
        projectResults.thetaModulation.phasedistro_lowRes_doubled(ii,:) = ...
            zscore(smooth([projectResults.thetaModulation.phasedistro_lowRes(ii,:) projectResults.thetaModulation.phasedistro_lowRes(ii,:)],5));
    end
       
    figure % 2.5 x 7
    subplot(4,1,[1 2])
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(responsive_cells,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(responsive_cells));
    ax = axis; xlim([0 4*pi]);
    x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave) + 1)/2) * ax(4)/3;
    hold on; plot(x_wave, y_wave,'color',[.3 .3 .3]); colormap jet
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel(name_cells,'Color',color_cells);
    subplot(4,1,3)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_pyr,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_pyr));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[]); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    subplot(4,1,4)
    imagesc_ranked([ps_lowRes; ps_lowRes+2*pi],[],...
        [projectResults.thetaModulation.phasedistro_lowRes_doubled(sessions_int,:)],[-3 3],...
        projectResults.thetaModulation.phasestats_m(sessions_int));
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('INT','Color',color_int);
    xlabel('Theta phase (rad)');
    colormap jet
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.thetaModulation.phasestats_m(sessions_pyr), projectResults.thetaModulation.phasestats_m(sessions_int),...
        projectResults.thetaModulation.phasestats_m(responsive_cells)},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType',...
        'roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Theta preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.thetaModulation.phasestats_r(sessions_pyr)), (projectResults.thetaModulation.phasestats_r(sessions_int)),...
        (projectResults.thetaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    
%     % 2. DIMENSIONALITY REDUCTION
%     Y = tsne([projectResults.cell_metrics.troughToPeak'...
%         projectResults.cell_metrics.acg_tau_rise'...
%         log10(projectResults.cell_metrics.firingRate_MA)'...
%         projectResults.cell_metrics.firingRateCV'...
%         projectResults.averageCCG.peakResponseZ...
%         projectResults.ripplesResponses.peakResponseZ_norm...
%         projectResults.ripplePhaseModulation.phasestats_m...
%         projectResults.ripplePhaseModulation.phasestats_r...
%         projectResults.thetaModulation.phasestats_m...
%         projectResults.thetaModulation.phasestats_r]);
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),30,color_cells_dark,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     
%     % plot features
%     figure
%     subplot(3,2,1)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7,log10(projectResults.cell_metrics.troughToPeak(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('troughToPeak','FontWeight','normal');
%     subplot(3,2,2)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.acg_tau_rise(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('acg_tau_rise','FontWeight','normal');
%     subplot(3,2,3)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.cell_metrics.firingRate_MA(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('firingRate_MA','FontWeight','normal');
%     xlabel('Dimension 2');
%     subplot(3,2,4)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, log10(projectResults.ripplesResponses.peakResponseZ_norm(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('ripple response','FontWeight','normal');
%     subplot(3,2,5)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.ripplePhaseModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple phase','FontWeight','normal'); xlabel('                                              Dimension 1');
%     subplot(3,2,6)
%     hold on
%     scatter(Y((targetSessCells),1),Y((targetSessCells),2),7, (projectResults.thetaModulation.phasestats_m(targetSessCells)),'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),15,color_cells_dark,'LineWidth',1);
%     colormap parula; title('Ripple response','FontWeight','normal');
    
    
    
    
%     % which are the other cells?
%     responsive_cells_allSessions = any(projectResults.optogeneticResponses.threeWaysTest'==1);
%     
%     targetSessCells_ivy = strcmpi(projectResults.geneticLine,'id2/nkx2.1/ai80');
%     targetSessCells_pv = strcmpi(projectResults.geneticLine,'pv/ai32');
%     targetSessCells_sst = strcmpi(projectResults.geneticLine,'sst/ai32');
%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),50,color_cells_dark,'LineWidth',2);
%     
%     scatter(Y(responsive_cells_allSessions(targetSessCells_ivy),1),Y(responsive_cells_allSessions(targetSessCells_ivy),2),50,color_id2nkx,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
%     scatter(Y(responsive_cells_allSessions(targetSessCells_pv),1),Y(responsive_cells_allSessions(targetSessCells_pv),2),50,color_pv,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
% %     scatter(Y(responsive_cells_allSessions(targetSessCells_sst),1),Y(responsive_cells_allSessions(targetSessCells_sst),2),50,color_sst,'LineWidth',2);
% %     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
end

% CHAPTER 5: CELLS COMPARISON
for z = 1
    % LIGHT RESPONSES
    [gs] = groupStats({projectResults.optogeneticResponses.rateDuringPulse(allResponsive_cells.id2,post_opt), projectResults.optogeneticResponses.rateDuringPulse(allResponsive_cells.pv,post_opt),...
        projectResults.optogeneticResponses.rateDuringPulse(allResponsive_cells.sst,post_opt),projectResults.optogeneticResponses.rateDuringPulse(allResponsive_cells.camk2,post_opt)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','symRoundPlot','labelSummary',false);
    ylabel('Light responses (Hz)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
    % SPIKE WIDHT
    cell_metrics = projectResults.cell_metrics;
    waveforms_timestmaps = cell_metrics.waveforms.time{1};
    
    figure
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,allResponsive_cells.id2),'color',color_id2dlx);
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,allResponsive_cells.pv),'color',color_pv);
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,allResponsive_cells.sst),'color',color_sst);
    plotFill(waveforms_timestmaps, cell_metrics.waveforms.filt_zscored(:,allResponsive_cells.camk2),'color',color_camk2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(allResponsive_cells.id2), projectResults.cell_metrics.troughToPeak(allResponsive_cells.pv),...
        projectResults.cell_metrics.troughToPeak(allResponsive_cells.sst),projectResults.cell_metrics.troughToPeak(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false);
    ylabel('Spike width (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
    % ACG
    acg_timestmaps = projectResults.acgPeak.acg_time(1,:);
    figure
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(allResponsive_cells.pv,:),'color',color_pv);
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(allResponsive_cells.sst,:),'color',color_sst);
    plotFill(acg_timestmaps, projectResults.acgPeak.acg_smoothed_norm(allResponsive_cells.camk2,:),'color',color_camk2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    [gs] = groupStats({projectResults.cell_metrics.troughToPeak(allResponsive_cells.id2), projectResults.cell_metrics.troughToPeak(allResponsive_cells.pv),...
        projectResults.cell_metrics.troughToPeak(allResponsive_cells.sst),projectResults.cell_metrics.troughToPeak(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false);
    ylabel('Spike width (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
    [gs] = groupStats({projectResults.acgPeak.peakTime(allResponsive_cells.id2), projectResults.acgPeak.peakTime(allResponsive_cells.pv),...
        projectResults.acgPeak.peakTime(allResponsive_cells.sst),projectResults.acgPeak.peakTime(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false);
    ylim(log10([0.001 1])); LogScale('y',10); ylim(log10([0.0015 1]));
    ylabel('ACG peak (s)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
    % QUIET VS RUN, NREM vs REM
    figure
    subplot(1,2,1)
    [gs] = groupStats({cell_metrics.run_quiet(allResponsive_cells.id2), cell_metrics.run_quiet(allResponsive_cells.pv),...
        cell_metrics.run_quiet(allResponsive_cells.sst),cell_metrics.run_quiet(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    hold on
    plot([0 5],[0 0],'color',[.7 .7 .7]); ylabel('NREM <----> REM'); ylim([-.5 .5]);
    ylabel('Spike width (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    ylabel('Quiet <----> Run'); ylim([-1 1]);
    
    subplot(1,2,2)
    [gs] = groupStats({cell_metrics.rem_nrem(allResponsive_cells.id2), cell_metrics.rem_nrem(allResponsive_cells.pv),...
        cell_metrics.rem_nrem(allResponsive_cells.sst),cell_metrics.rem_nrem(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    hold on
    plot([0 5],[0 0],'color',[.7 .7 .7]); ylabel('NREM <----> REM'); ylim([-.5 .5]);
    ylabel('Spike width (ms)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    ylabel('NREM <----> REM'); ylim([-1 1]);
    
    
   % CCG
    ts_CCG = projectResults.averageCCG.timestamps;
    figure
    plotFill(ts_CCG, projectResults.averageCCG.ZmeanCCGSmooth(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill(ts_CCG, projectResults.averageCCG.ZmeanCCGSmooth(allResponsive_cells.pv,:),'color',color_pv);
    plotFill(ts_CCG, projectResults.averageCCG.ZmeanCCGSmooth(allResponsive_cells.sst,:),'color',color_sst);
    plotFill(ts_CCG, projectResults.averageCCG.ZmeanCCGSmooth(allResponsive_cells.camk2,:),'color',color_camk2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    [gs] = groupStats({projectResults.averageCCG.peakResponseZ(allResponsive_cells.id2), projectResults.averageCCG.peakResponseZ(allResponsive_cells.pv),...
        projectResults.averageCCG.peakResponseZ(allResponsive_cells.sst),projectResults.averageCCG.peakResponseZ(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false);
    ylabel('AvgCCG (Z)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
   
   % THETA
    ps = projectResults.thetaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    figure
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.thetaModulation.phasedistro_lowRes_doubled(allResponsive_cells.camk2,:),'color',color_camk2);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.thetaModulation.phasedistro_lowRes_doubled(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.thetaModulation.phasedistro_lowRes_doubled(allResponsive_cells.pv,:),'color',color_pv);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.thetaModulation.phasedistro_lowRes_doubled(allResponsive_cells.sst,:),'color',color_sst);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]); ylabel('PYR','Color',color_pyr);
    ylabel('Rate (SD)');
    xlabel('Theta phase (rad)');
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.thetaModulation.phasestats_m(allResponsive_cells.id2), projectResults.thetaModulation.phasestats_m(allResponsive_cells.pv),...
        projectResults.thetaModulation.phasestats_m(allResponsive_cells.sst),projectResults.thetaModulation.phasestats_m(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Theta preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({projectResults.thetaModulation.phasestats_r(allResponsive_cells.id2), projectResults.thetaModulation.phasestats_r(allResponsive_cells.pv),...
        projectResults.thetaModulation.phasestats_r(allResponsive_cells.sst),projectResults.thetaModulation.phasestats_r(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
   % GAMMA
    ps = projectResults.lGammaModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    figure
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.lGammaModulation.phasedistro_lowRes_doubled(allResponsive_cells.camk2,:),'color',color_camk2);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.lGammaModulation.phasedistro_lowRes_doubled(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.lGammaModulation.phasedistro_lowRes_doubled(allResponsive_cells.pv,:),'color',color_pv);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.lGammaModulation.phasedistro_lowRes_doubled(allResponsive_cells.sst,:),'color',color_sst);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]);
    ylabel('Rate (SD)');
    xlabel('Gamma phase (rad)');
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.lGammaModulation.phasestats_m(allResponsive_cells.id2), projectResults.lGammaModulation.phasestats_m(allResponsive_cells.pv),...
        projectResults.lGammaModulation.phasestats_m(allResponsive_cells.sst),projectResults.lGammaModulation.phasestats_m(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Gamma phase preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({projectResults.lGammaModulation.phasestats_r(allResponsive_cells.id2), projectResults.lGammaModulation.phasestats_r(allResponsive_cells.pv),...
        projectResults.lGammaModulation.phasestats_r(allResponsive_cells.sst),projectResults.lGammaModulation.phasestats_r(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
   
   
   % RIPPLE
    ts_ripples = projectResults.ripplesResponses.timestamps;
    figure
    plotFill(ts_ripples, projectResults.ripplesResponses.responseZ(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill(ts_ripples, projectResults.ripplesResponses.responseZ(allResponsive_cells.pv,:),'color',color_pv);
    plotFill(ts_ripples, projectResults.ripplesResponses.responseZ(allResponsive_cells.sst,:),'color',color_sst);
    plotFill(ts_ripples, projectResults.ripplesResponses.responseZ(allResponsive_cells.camk2,:),'color',color_camk2);
    axis tight; xlabel('Time (ms)'); ylabel('Ripple response (SD)'); set(gca, 'TickDir', 'out'); xlim([-.2 .2])
    
    [gs] = groupStats({projectResults.ripplesResponses.peakResponseZ(allResponsive_cells.id2), projectResults.ripplesResponses.peakResponseZ(allResponsive_cells.pv),...
        projectResults.ripplesResponses.peakResponseZ(allResponsive_cells.sst),projectResults.ripplesResponses.peakResponseZ(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false);
    ylabel('AvgCCG (Z)');    
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ID2+/DLX+','PV+','SST+','CaMKII'},'XTickLabelRotation',45);
    
    
    % RIPPLE PHASE
    ps = projectResults.ripplePhaseModulation.phasebins(1,:);
    ps_bins = discretize(ps,36);
    ps_lowRes = wrapTo2Pi(accumarray(ps_bins',ps,[],@mean));
    
    figure
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(allResponsive_cells.camk2,:),'color',color_camk2);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(allResponsive_cells.id2,:),'color',color_id2dlx);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(allResponsive_cells.pv,:),'color',color_pv);
    plotFill([ps_lowRes; ps_lowRes+2*pi], projectResults.ripplePhaseModulation.phasedistro_lowRes_doubled(allResponsive_cells.sst,:),'color',color_sst);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    ax = axis; xlim([0 4*pi]);
    set(gca,'XTick',[],'XTick',[0:pi:4*pi],'XTickLabel',{'0','π','2π','3π','4π'}); xlim([0 4*pi]);
    ylabel('Rate (SD)');
    xlabel('Ripple phase (rad)');
    
    figure
    subplot(1,2,1)
    [gs] = groupStats({projectResults.ripplePhaseModulation.phasestats_m(allResponsive_cells.id2), projectResults.ripplePhaseModulation.phasestats_m(allResponsive_cells.pv),...
        projectResults.ripplePhaseModulation.phasestats_m(allResponsive_cells.sst),projectResults.ripplePhaseModulation.phasestats_m(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'orientation','horizontal','inAxis',true);
    hold on
    ax = axis; ylim([0 2*pi]);
    x_wave = deg2rad([0:0.5:360]); y_wave =-(((cos(x_wave) + 1)) * ax(2)/3) + 3.5;
    plot(y_wave, x_wave,'color',[.3 .3 .3]);
    set(gca,'XTick',[],'YTick',[0:pi:4*pi],'YTickLabel',{'0','π','2π','3π','4π'},'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells});
    ylabel('Ripple phase preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({projectResults.ripplePhaseModulation.phasestats_r(allResponsive_cells.id2), projectResults.ripplePhaseModulation.phasestats_r(allResponsive_cells.pv),...
        projectResults.ripplePhaseModulation.phasestats_r(allResponsive_cells.sst),projectResults.ripplePhaseModulation.phasestats_r(allResponsive_cells.camk2)},...
        [],'color',[color_id2dlx; color_pv; color_sst; color_camk2],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Mean Vector Length');    
    set(gca,'XTick',[1 2 3],'XTickLabel',{'PYR','INT',name_cells},'XTickLabelRotation',45);
    
    
    % tSNE
    allResponsive_cells.all = [find(allResponsive_cells.id2) find(allResponsive_cells.pv) find(allResponsive_cells.sst) find(allResponsive_cells.camk2)]';
    allResponsive_cells.all_group = [ones(length(find(allResponsive_cells.id2)),1); 2*ones(length(find(allResponsive_cells.pv)),1); 3*ones(length(find(allResponsive_cells.sst)),1); 4*ones(length(find(allResponsive_cells.camk2)),1)];

    projectResults.cell_metrics.troughToPeak
    
    
    Y = tsne([projectResults.cell_metrics.troughToPeak(allResponsive_cells.all)' ...
        projectResults.acgPeak.peakTime(allResponsive_cells.all)' ...
        projectResults.thetaModulation.phasestats_m(allResponsive_cells.all) ...
        projectResults.thetaModulation.phasestats_r(allResponsive_cells.all)]);
    
    
    figure;
    hold on
    scatter(Y(allResponsive_cells.all_group==1,1),Y(allResponsive_cells.all_group==1,2),30,color_id2dlx,'filled');
    scatter(Y(allResponsive_cells.all_group==2,1),Y(allResponsive_cells.all_group==2,2),30,color_pv,'filled');
    scatter(Y(allResponsive_cells.all_group==3,1),Y(allResponsive_cells.all_group==3,2),30,color_sst,'filled');
    scatter(Y(allResponsive_cells.all_group==4,1),Y(allResponsive_cells.all_group==4,2),30,color_camk2,'filled');
    ylabel('Dim2'); xlabel('Dim1');
    
    
    figure;
    hold on
    scatter(Y(allResponsive_cells.all_group==1,1),Y(allResponsive_cells.all_group==1,2),30,color_id2dlx,'filled');
    scatter(Y(allResponsive_cells.all_group==2,1),Y(allResponsive_cells.all_group==2,2),30,color_pv,'filled');
    scatter(Y(allResponsive_cells.all_group==3,1),Y(allResponsive_cells.all_group==3,2),30,color_sst,'filled');
    scatter(Y(allResponsive_cells.all_group==4,1),Y(allResponsive_cells.all_group==4,2),30,color_camk2,'filled');
    ylabel('Dim2'); xlabel('Dim1');

%     
%     
%     figure
%     hold on
%     scatter(Y(is_pyr(targetSessCells),1),Y(is_pyr(targetSessCells),2),10,color_pyr_light,'filled');
%     scatter(Y(is_int(targetSessCells),1),Y(is_int(targetSessCells),2),10,color_int_light,'filled');
%     scatter(Y(responsive_cells(targetSessCells),1),Y(responsive_cells(targetSessCells),2),30,color_cells_dark,'LineWidth',2);
%     xlabel('Dimension 1');  ylabel('Dimension 2'); title('tSNE','FontWeight','normal');
    
end



