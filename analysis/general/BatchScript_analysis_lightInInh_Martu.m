%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear all; close all
targetProject= 'LightInInh';
targetBehavior = 'linear maze';

HCB_directory = what('HippoCookBook'); 


sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions


for ii = 183: length(sessionsTable.SessionName)
     %% Analysis general all over Camkii/32 animal
    if contains(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)

        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}]));

        % try
            %%% your code goes here...
            clear uLEDResponses_interval
            delete(gcp('nocreate'))
            spikes = loadSpikes;
            spikes_times = spikes.times;
            monosyn_inh_win = [0.001 0.021]; % 01to21
            monosyn_control_win = [-0.041 -0.021];

            for mm = 1:spikes.numcells 
                disp(mm);
                uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
                    'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false, 'boostraping_type','pulses');
            end

            collision_metrics_1_21 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','1msTo21ms','saveMat',true,'update_cell_metrics',true,'save_as','lightSpikeCollisions','rate_change_threshold',3);

            clear uLEDResponses_interval
            clear uLEDResponses_control
            close all;
       
    end
end


%% Analysis for pre and post synaptic changes


for ii = 86:length(sessionsTable.SessionName)

    if contains(sessionsTable.Project{ii}, targetProject) && contains(sessionsTable.Behavior{ii},targetBehavior) || strcmpi('all', targetProject) 
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}])); 

        session = loadSession;
        delete(gcp('nocreate'))
        spikes = loadSpikes;
        spikes_times = spikes.times;
        monosyn_inh_win = [0.001 0.021]; % 01to21

        for kk= 1:length(session.epochs) 
            if strcmp(session.epochs{kk}.behavioralParadigm,'Maze')
                pre_end = session.epochs{kk}.startTime;
                post_start = session.epochs{kk}.stopTime;
                pre_maze = [0 pre_end];
                post_maze = [post_start Inf];
            end
        end
 
        for mm = 1:spikes.numcells
            disp(mm);

            uLEDResponses_interval_pre{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
                'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false,'restrict_to',pre_maze,'minNumberOfPulses',20);

            uLEDResponses_interval_post{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
                'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false,'restrict_to',post_maze, 'minNumberOfPulses',20);
        end

        collision_metrics_1_21_pre = get_light_spike_CollisionMetrics(uLEDResponses_interval_pre,'label','1msTo21ms_pre','saveMat',true,'update_cell_metrics',false,'save_as','lightSpikeCollisions_pre','rate_change_threshold',3);

        collision_metrics_1_21_post = get_light_spike_CollisionMetrics(uLEDResponses_interval_post,'label','1msTo21ms_post','saveMat',true,'update_cell_metrics',false,'save_as','lightSpikeCollisions_post','rate_change_threshold',3);
        %%%

        close all;
        clear uLEDResponses_interval_pre
        clear uLEDResponses_interval_post


    end
end

%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

% 
% session = loadSession;
%  for ii= 1:length(session.epochs) 
%     if strcmp(session.epochs{ii}.behavioralParadigm,'Maze')
%         pre_end = session.epochs{ii}.startTime;
%         post_start = session.epochs{ii}.stopTime;
%         pre_maze = [0 pre_end];
%         post_maze = [post_start Inf];
%     end
%  end
%         [status] = InIntervals(timestamps,post_maze);
%         sum(status)

    for ii = 1: length(sessionsTable.SessionName)
      
        if contains(sessionsTable.Project{ii}, targetProject) && contains(sessionsTable.Behavior{ii},targetBehavior) || strcmpi('all', targetProject) 
    
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}])); 
    
        session = loadSession;
        delete(gcp('nocreate'))
        spikes = loadSpikes;
        spikes_times = spikes.times;
        monosyn_inh_win = [0.001 0.021]; % 01to21
    
        for kk= 1:length(session.epochs) 
            if strcmp(session.epochs{kk}.behavioralParadigm,'Maze')
                pre_end = session.epochs{kk}.startTime;
                post_start = session.epochs{kk}.stopTime;
                pre_maze = [0 pre_end];
                post_maze = [post_start Inf];
            end
        end
    
        psthRipples_before = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to',pre_maze,'save_as','psthRipples_before');
        psthRipples_after = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to',post_maze,'save_as','psthRipples_after');
        psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);

        close all

        cell_metrics = ProcessCellMetrics('session', session,'getWaveformsFromDat',false,'forceReload',true);
        cell_metrics_before = ProcessCellMetrics('session', session,'getWaveformsFromDat',false,'forceReload',true,'restrictToIntervals',pre_maze,'excludeMetrics',{'deepSuperficial','monoSynaptic_connections'},'saveAs','cell_metrics_before');
        cell_metrics_after = ProcessCellMetrics('session', session,'getWaveformsFromDat',false,'forceReload',true,'restrictToIntervals',post_maze,'excludeMetrics',{'deepSuperficial','monoSynaptic_connections'},'saveAs','cell_metrics_after'); % after CellExplorar
       
        close all
    
        end
    end 
