

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project_script_interneuronsLibrary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MV, PA, RGH 2022
% figures in figures_script_interneuronsLibrary

% CHAPTER_0: GET DATA
for z = 1
    clear; close all
    analysis_project_path = adapt_filesep([dropbox_path filesep 'ProjectsOnLine\interneuronsLibrary\data']);
    [projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary','analysis_project_path', analysis_project_path,'loadLast',true);
    
    % general
    inHippocampus = {'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1' 'CA3' 'DG'}; % only using hippocampus data... :)
    
    is_pyr = strcmpi(projectResults.cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    is_int = (strcmpi(projectResults.cell_metrics.putativeCellType,'Narrow Interneuron')...
        | strcmpi(projectResults.cell_metrics.putativeCellType,'Wide Interneuron')) & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    color_id2dlx = [255 180 40]/255;
    color_id2dlx_dark = [155 90 20]/255;
    
    color_id2nkx = [200 200 40]/255;
    color_id2nkx_dark = [155 155 20]/255;
    
    color_sst = [.3 .9 .3];
    color_sst_dark = [.2 .6 .2];
    
    color_pv = [.3 .3 .9];
    color_pv_dark = [.2 .2 .6];
    
    color_camk2 = [.9 .5 .3];
    color_camk2_dark = [.6 .3 .1];
    
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
    projectResults.ripplesResponses.peakResponseZ_norm = ZeroToOne(projectResults.ripplesResponses.peakResponseZ);
end

% CHAPTER_1: ID2 Interneurons (Neurogliaform and CCK basket
for z = 1
    % ID2/Dlx cells features
    targetSessCells = strcmpi(projectResults.geneticLine,'id2/dlx/ai80') & ismember(projectResults.cell_metrics.brainRegion,inHippocampus);
    allTargetSessCells.id2 = targetSessCells;
    
    name_cells = 'ID2+/DLX+';
    color_cells = color_id2dlx;
    color_cells_dark = color_id2dlx_dark;
    
    % 1.1 Light responses
    ts = projectResults.optogeneticResponses.timestamps;
    responsive_cells = any(projectResults.optogeneticResponses.threeWaysTest'==1) & targetSessCells & (projectResults.optogeneticResponses.checkedCells==1)';
    responsive_cells([137 945]) = 0;
    
    allResponsive_cells.id2 = responsive_cells;
    [~,optimal_pulse] = max(projectResults.optogeneticResponses.rateDuringPulse');
    sessions_pyr = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_pyr;
    sessions_int = ~any(projectResults.optogeneticResponses.threeWaysTest') & targetSessCells & is_int;
    all_pyr.id2 = sessions_pyr;
    all_int.id2 = sessions_int;
    
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
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[-0.1],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(sessions_int)), log10(cell_metrics.firingRate_WAKEnontheta(sessions_int)),...
        log10(cell_metrics.firingRate_NREMstate(sessions_int)), log10(cell_metrics.firingRate_REMstate(sessions_int))},[],'color',[color_int],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0],'sigStar',false);
    [gs] = groupStats({log10(cell_metrics.firingRate_WAKEtheta(responsive_cells)), log10(cell_metrics.firingRate_WAKEnontheta(responsive_cells)),...
        log10(cell_metrics.firingRate_NREMstate(responsive_cells)), log10(cell_metrics.firingRate_REMstate(responsive_cells))},[],'color',[color_cells],...
        'plotData',true,'plotType','meanBallCI95','labelSummary',false,'inAxis',true,'posOffset',[0.1],'sigStar',false);
    ylim(log10([0.05 40])); LogScale('y',10);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Run','Quiet','NREM','REM'},'XTickLabelRotation',45);
    ylim(log10([0.3 5])); LogScale('y',10);
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
    
    % 1.5 THETA
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
    ylabel('GammaM preference (deg)');
    subplot(1,2,2)
    [gs] = groupStats({(projectResults.lGammaModulation.phasestats_r(sessions_pyr)), (projectResults.lGammaModulation.phasestats_r(sessions_int)),...
        (projectResults.lGammaModulation.phasestats_r(responsive_cells))},[],'color',[color_pyr; color_int; color_cells],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
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



