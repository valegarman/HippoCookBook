%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'LightInInh';

HCB_directory = what('HippoCookBook'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    if contains(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{ii}]));
        try
        
            %%% your code goes here...
            clear uLEDResponses_interval
            delete(gcp('nocreate'))
            spikes = loadSpikes;
            spikes_times = spikes.times;
            monosyn_inh_win = [0.01 0.021]; % 01to21
            parfor mm = 1:spikes.numcells
                disp(mm);
                uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
                    'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
            end
            collision_metrics_01_21 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','01msTo21ms','saveMat',false,'update_cell_metrics',false);
            save('uLEDResponses_interval_01ms_21ms.mat','uLEDResponses_interval','collision_metrics_01_21');
            clear uLEDResponses_interval

            load('uLEDResponses_interval_01ms_21ms.mat');
            collision_metrics_01_21 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','01msTo21ms','rate_change_threshold',3);

            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end

%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions