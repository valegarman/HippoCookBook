%% BatchScript_analysis_ripplesPerSubsession

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'POTometry';
force_fiber_pre_post = true;

% Indexes for cell location
inHippocampus = {'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1' 'CA3' 'DG' 'CA3sp' 'CA3sr'}; % only using hippocampus data... :);
inCA1 = {'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1'};
inForebrain = {'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1' 'CA3' 'DG','CA3sp' 'CA3sr' 'PTLp' 'PTLp2_3' 'PTLp4' 'PTLp5' 'PTLp6', 'VISp1', 'VISp2_3', 'VISp4', 'VISp5', 'VISp6'};
inCortex = {'PTLp' 'PTLp2_3' 'PTLp4' 'PTLp5' 'PTLp6', 'VISp1', 'VISp2_3', 'VISp4', 'VISp5', 'VISp6'};

for ii = 1:length(sessionsTable.SessionName)
    if contains(sessionsTable.Project(ii), targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}]);
        basepath = pwd;

        % 
        clear basepath session basename cell_metrics spikes 
        basepath = pwd;
        session = loadSession;
        basename = basenameFromBasepath(basepath);
        spikes = loadSpikes();
        % Load cell_metrics
        file = dir(['*cell_metrics.cellinfo.mat']);load(file.name);

        % Indexes for predicted cell type
        % is_pyr = strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(cell_metrics.brainRegion,inHippocampus);
        is_pyr = strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell');
        is_pyr_CA1 = strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(cell_metrics.brainRegion,inCA1);
        
        % is_Camk2 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_Camk2 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')';
        is_Camk2_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_Camk2_deep = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_Camk2_deep = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')';
        is_Camk2_deep_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_Camk2_sup = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_Camk2_sup = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')';
        is_Camk2_sup_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_pv = strcmpi(cell_metrics.ground_truth_classification.cell_types,'PV+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_pv = strcmpi(cell_metrics.ground_truth_classification.cell_types,'PV+')';
        is_pv_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'PV+') & ismember(cell_metrics.brainRegion,inCA1);

        % is_sst = strcmpi(cell_metrics.ground_truth_classification.cell_types,'SST+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_sst = strcmpi(cell_metrics.ground_truth_classification.cell_types,'SST+')';
        is_sst_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'SST+')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_id2 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'Id2+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_id2 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'Id2+')';
        is_id2_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'Id2+') & ismember(cell_metrics.brainRegion,inCA1);

        % is_vip = strcmpi(cell_metrics.ground_truth_classification.cell_types,'VIP+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_vip = strcmpi(cell_metrics.ground_truth_classification.cell_types,'VIP+')';
        is_vip_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_types,'VIP+')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_no_sncg = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_no_sncg = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')';
        is_no_sncg_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')' & ismember(cell_metrics.brainRegion,inCA1);

        % is_sncg = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')' & ismember(cell_metrics.brainRegion,inHippocampus);
        is_sncg = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')';
        is_sncg_CA1 = strcmpi(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')' & ismember(cell_metrics.brainRegion,inCA1);

        % if isempty(dir('*fiber_psth_ripples_PreSleep2.mat'))
            % 
            % ts_PreSleep2 = [];
            % ts_PostSleep2 = [];
            % ripples_fiber_pre = [];
            % ripples_fiber_post = [];
            % 
            % for jj = 1:length(session.epochs)
            %     if strcmpi(session.epochs{jj}.behavioralParadigm,'PreSleep2') && strcmpi(session.epochs{jj}.manipulation,'Fiber')
            %         ts_PreSleep2 = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            %         path_PreSleep2 = session.epochs{jj}.name;
            %     elseif strcmpi(session.epochs{jj}.behavioralParadigm,'PreSleep2') && ~strcmpi(session.epochs{jj}.manipulation,'Fiber')
            %         error('Problems defining fiber epochs...');
            %     end
            % end
            % 
            % for jj = 1:length(session.epochs)
            %     if strcmpi(session.epochs{jj}.behavioralParadigm,'PostSleep2') && strcmpi(session.epochs{jj}.manipulation,'Fiber')
            %         ts_PostSleep2 = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            %         path_PostSleep2 = session.epochs{jj}.name;
            %     elseif strcmpi(session.epochs{jj}.behavioralParadigm,'PostSleep2') && ~strcmpi(session.epochs{jj}.manipulation,'Fiber')
            %         error('Problems defining fiber epochs...');  
            %     end
            % end

            % Compute fiber for pre-post maze epochs
            % if ~isempty(ts_PreSleep2)
            %     cd(path_PreSleep2)
            %     peri_spike_trace = peri_spike_trace_average(spikes);
            %     cd(basepath);
            %     close all;
            %     save([session.general.name,'.peri_spike_trace_PreSleep2.mat'],'peri_spike_trace');
            % 
            %     % plot
            %     if length(find(is_Camk2)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_Camk2,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_Camk2,:));
            %     end
            %     if length(find(is_pv)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_pv,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_pv,:));
            %     end
            % 
            %     if length(find(is_sst)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_sst,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_sst,:));
            %     end
            % 
            %     if length(find(is_id2)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_id2,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_id2,:));
            %     end
            % 
            %     if length(find(is_vip)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_vip,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_vip,:));
            %     end
            % 
            % 
            % end
            % 
            % peri_spike_trace = [];
            % if ~isempty(ts_PostSleep2)
            %     cd(path_PostSleep2);
            %     peri_spike_trace = peri_spike_trace_average(spikes);     
            %     cd(basepath);
            %     close all;
            %     save([session.general.name,'.peri_spike_trace_PostSleep2.mat'],'peri_spike_trace');
            % 
            %     close all;
            %     % plot
            %     if length(find(is_Camk2)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_Camk2,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_Camk2,:));
            %     end
            %     if length(find(is_pv)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_pv,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_pv,:));
            %     end
            % 
            %     if length(find(is_sst)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_sst,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_sst,:));
            %     end
            % 
            %     if length(find(is_id2)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_id2,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_id2,:));
            %     end
            % 
            %     if length(find(is_vip)) > 0
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_red_smoothed(is_vip,:));
            % 
            %         figure;
            %         plotFill(peri_spike_trace.lag_time(1,:),peri_spike_trace.corr_green_smoothed(is_vip,:));
            %     end
            % 
            % end
        % end

        close all;

        % 
        % Cell classifier
        % computeRippleReversal;
        % cellTypeClassifier;

         try
            ripples_fiber = fiberPhotometryModulation_temp([],'eventType','ripples','reload_fiber',true);
        catch
            warning('No fiber recording in this session...');
        end

        if isempty(dir('*fiber_psth_ripples_PreSleep2.mat')) || force_fiber_pre_post

            % 

            ts_PreSleep2 = [];
            ts_PostSleep2 = [];
            ripples_fiber_pre = [];
            ripples_fiber_post = [];

            for jj = 1:length(session.epochs)
                if strcmpi(session.epochs{jj}.behavioralParadigm,'PreSleep2') && strcmpi(session.epochs{jj}.manipulation,'Fiber')
                    ts_PreSleep2 = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                elseif strcmpi(session.epochs{jj}.behavioralParadigm,'PreSleep2') && ~strcmpi(session.epochs{jj}.manipulation,'Fiber')
                    error('Problems defining fiber epochs...');
                end
            end

            for jj = 1:length(session.epochs)
                if strcmpi(session.epochs{jj}.behavioralParadigm,'PostSleep2') && strcmpi(session.epochs{jj}.manipulation,'Fiber')
                    ts_PostSleep2 = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                elseif strcmpi(session.epochs{jj}.behavioralParadigm,'PostSleep2') && ~strcmpi(session.epochs{jj}.manipulation,'Fiber')
                    error('Problems defining fiber epochs...');
                end
            end

            % Compute fiber for pre-post maze epochs
            if ~isempty(ts_PreSleep2)
                ripples_fiber_pre = fiberPhotometryModulation_temp([],'restrictIntervals',ts_PreSleep2,'eventType','ripples','saveAs','PreSleep2','savePlotAs','PreSleep2');
            end
            if ~isempty(ts_PostSleep2)
                ripples_fiber_post = fiberPhotometryModulation_temp([],'restrictIntervals',ts_PostSleep2,'eventType','ripples','saveAs','PostSleep2','savePlotAs','PostSleep2');
            end

            try
                figure;
                plotFill([1:size(ripples_fiber_pre.red_normalized.responsecurveZSmooth,2)],ripples_fiber_pre.red_normalized.responsecurveZSmooth,'color',[0 0 0]);
                hold on;
                plotFill([1:size(ripples_fiber_post.red_normalized.responsecurveZSmooth,2)],ripples_fiber_post.red_normalized.responsecurveZSmooth,'color',[1 0 0]);
                saveas(gca,['SummaryFigures\fiber_red_pre_post.png']);

                figure;
                plotFill([1:size(ripples_fiber_pre.green_normalized.responsecurveZSmooth,2)],ripples_fiber_pre.green_normalized.responsecurveZSmooth,'color',[0 0 0]);
                hold on;
                plotFill([1:size(ripples_fiber_post.green_normalized.responsecurveZSmooth,2)],ripples_fiber_post.green_normalized.responsecurveZSmooth,'color',[0 1 0]);
                saveas(gca,['SummaryFigures\fiber_green_pre_post.png']);

                figure;
                plotFill([1:size(ripples_fiber_pre.red_PP.responsecurveZSmooth,2)],ripples_fiber_pre.red_PP.responsecurveZSmooth,'color',[0 0 0]);
                hold on;
                plotFill([1:size(ripples_fiber_post.red_PP.responsecurveZSmooth,2)],ripples_fiber_post.red_PP.responsecurveZSmooth,'color',[1 0 0]);
                saveas(gca,['SummaryFigures\fiber_red_PP_pre_post.png']);

                figure;
                plotFill([1:size(ripples_fiber_pre.green_PP.responsecurveZSmooth,2)],ripples_fiber_pre.green_PP.responsecurveZSmooth,'color',[0 0 0]);
                hold on;
                plotFill([1:size(ripples_fiber_post.green_PP.responsecurveZSmooth,2)],ripples_fiber_post.green_PP.responsecurveZSmooth,'color',[0 1 0]);
                saveas(gca,['SummaryFigures\fiber_green_pre_post.png']);

            catch
            end
        end

        close all;

        
    end
    close all;
end