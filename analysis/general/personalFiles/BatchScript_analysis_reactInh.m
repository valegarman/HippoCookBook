%% BatchScript_analysis_reactInh
% place your code to run an analysis across all sessions for a given
% project

%% analysis for figure 2

clear; close all
targetProject= 'All';
list_of_sessions = {'fCamk1_200827_sess9', 'fCamk1_200901_sess12', 'fCamk1_200902_sess13',...
    'fCamk1_200904_sess15','fCamk1_200908_sess16','fCamk1_200909_sess17','fCamk1_200910_sess18',...
    'fCamk1_200911_sess19','fCamk3_201028_sess10_cleanned','fCamk3_201029_sess11_cleanned',...
    'fCamk3_201030_sess12','fCamk3_201102_sess13','fCamk3_201103_sess14','fCamk3_201111_sess20',...
    'fCamk3_201113_sess22','fCamk3_201105_sess16','fCamk3_201106_sess17','fCamk3_201110_sess19',...
    'fCamk3_201109_sess18'};

HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            clear uLEDResponses_interval
            delete(gcp('nocreate'))
            ripples = rippleMasterDetector;
            padding = .05;
            uLEDResponses_ripples = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples',...
                    'minNumberOfPulses',5);

            session = loadSession;
            pre_win = [];
            post_win = [];
            for jj = 1:length(session.epochs)
                if strcmpi(session.epochs{jj}.behavioralParadigm,'Maze')
                    pre_win = [0 session.epochs{jj}.startTime];
                    post_win = [session.epochs{jj}.stopTime NaN];
                elseif strcmpi(session.epochs{jj}.behavioralParadigm,'BaselinePost')
                    post_win(2) = [session.epochs{jj}.stopTime];
                end
            end
                    
            uLEDResponses_ripples_pre = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples_pre',...
                    'minNumberOfPulses',5,'restrict_to',pre_win);

            uLEDResponses_ripples_post = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples_post',...
                    'minNumberOfPulses',5,'restrict_to',post_win);

            %%%
            groupStats({uLEDResponses_ripples_pre.out_interval.maxRespLED.rate(uLEDResponses_ripples_pre.stronglyDrivenCells) - uLEDResponses_ripples_pre.out_interval.maxRespLED.rateBeforePulse(uLEDResponses_ripples_pre.stronglyDrivenCells), ...
                uLEDResponses_ripples_pre.in_interval.maxRespLED.rate(uLEDResponses_ripples_pre.stronglyDrivenCells) - uLEDResponses_ripples_pre.in_interval.maxRespLED.rateBeforePulse(uLEDResponses_ripples_pre.stronglyDrivenCells)})
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
end


%% analysis for figure 4
clear; close all
targetProject= 'All';
list_of_sessions = {'fcamk10_220921_sess9', 'fcamk10_220922_sess10', 'fcamk10_220927_sess13', 'fcamk10_220928_sess14', 'fcamk10_220929_sess15','fcamk10_221004_sess18'};

HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            targetFile = dir('*spikeTriggeredPulses.cellinfo.mat'); load(targetFile.name);
            spikes = loadSpikes;
            spikeTriggeredPulses.units_in_stim_shank = spikes.shankID == spikeTriggeredPulses.directStimShank;
            spikeTriggeredPulses.units_in_source_shank = spikes.shankID == spikeTriggeredPulses.controlShank;
            spikeTriggeredPulses.units_in_delayed_stim_shank = spikes.shankID == spikeTriggeredPulses.delayedStimShank(end);
            spikeTriggeredPulses.triggered_unit = zeros(size(spikes.shankID));
            spikeTriggeredPulses.triggered_unit(spikeTriggeredPulses.unitsTriggeringPulse(1)) = 1;
            save(targetFile.name,'spikeTriggeredPulses');

            uledPulses = getuLEDPulses;
            uledResponses = getuLEDResponse;
            stimulationEpoch = spikeTriggeredPulses.stimulationInterval;

            % explained variance of cells in stimulated shanks
            spikes = loadSpikes;
            spikes.times(~(spikes.shankID==spikeTriggeredPulses.controlShank) & ~(spikes.shankID==spikeTriggeredPulses.directStimShank)) = [];
            ev = explained_variance(spikes, [0 stimulationEpoch(1)],[stimulationEpoch],[stimulationEpoch(2) Inf],'save_as','explained_variance_stim');
            
            % explained variance of cells non stimulated shanks
            spikes = loadSpikes;
            spikes.times(~(spikes.shankID==spikeTriggeredPulses.controlShank) & ~(spikes.shankID==spikeTriggeredPulses.delayedStimShank(end))) = [];
            ev = explained_variance(spikes, [0 stimulationEpoch(1)],[stimulationEpoch],[stimulationEpoch(2) Inf],'save_as','explained_variance_delayed');

            % 
            [spikeCCGchange] = getSpikeCCGchange(spikeTriggeredPulses.unitsTriggeringPulse,[0 stimulationEpoch(1); stimulationEpoch(1) Inf; stimulationEpoch],'force',true);

            % uledresponses before and after
            padding = [-0.02 0.02];
            spikes = loadSpikes;

            uLEDResponses_spikeTriggere_post = getuLEDResponse_intervals([spikes.times{spikeTriggeredPulses.unitsTriggeringPulse(1)}+padding(1) spikes.times{spikeTriggeredPulses.unitsTriggeringPulse(1)}+padding(2)],...
                    'saveMat', true,'numRep',10,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_spikeTriggered',...
                    'minNumberOfPulses',5,'restrict_to',[stimulationEpoch(2) Inf],'doPlot',false);
            close all;
        catch
            warning('Analysis was not possible!');
        end
end

%% Analysis for co-activation
clear; close all
targetProject= 'All';
list_of_sessions = {'fCamk10_220913_sess3', 'fCamk10_220914_sess4', 'fCamk10_220915_sess5', 'fCamk10_220916_sess6', 'fCamk10_220930_sess16', 'fCamk10_221007_sess21'};
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{targetSessions(ii)}) filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            % identify times before, during, and after experiment
            getuLEDCoactivation('winCoactivation', 0.010);
            % 
            
        catch
            warning('Analysis was not possible!');
        end
end

% to take a baseline session before and after uled stim, I am going to take
% baselinePre, PostStim, BaselinePost
% 
list_of_sessions_control = {'fcamk1_200827_sess9', 'fcamk1_200901_sess12', 'fcamk1_200902_sess13', 'fcamk1_200904_sess15', 'fcamk1_200908_sess16', 'fcamk1_200909_sess17', 'fcamk1_200911_sess19', 'fcamk1_200910_sess18', ...
    'fcamk10_220920_sess8', 'fcamk3_201117_sess24', 'fcamk3_201030_sess12', 'fcamk3_201103_sess14', 'fcamk3_201111_sess20', 'fcamk3_201105_sess16','fcamk3_201109_sess18', 'fcamk3_201029_sess11_cleanned', 'fcamk3_201102_sess13', ...
    'fcamk3_201110_sess19', 'fcamk3_201113_sess22', 'fcamk3_201106_sess17', 'fcamk3_201116_sess23'};
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));
for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{targetSessions(ii)}) filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            % identify times before, during, and after experiment
            getuLEDCoactivation('winCoactivation', 0.010);
            % 
            
        catch
            warning('Analysis was not possible!');
        end
end



[projectResults, projectSessionResults] = ...
        loadProjectResults('list_of_sessions',list_of_sessions_control,...
        'analysis_project_path', [onedrive_path 'NeuralComputationLab\ActiveProjects\ReactInh\dataCoactivation'] ,'loadLast',true, 'save_as', 'reactInh_figureCoactivation_control');

% Control for coactivation
clear; close all
targetProject= 'All';
list_of_sessions = {'fCamk10_220913_sess3', 'fCamk10_220914_sess4', 'fCamk10_220915_sess5', 'fCamk10_220916_sess6', 'fCamk10_220930_sess16', 'fCamk10_221007_sess21'};


HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));


% 1. Arreglar el error
% 2. Loop por los pares
    list_of_results = {'optogeneticResponse','averageCCG','ripples_psth','theta_*.PhaseLockingData',...
            'thetaRun*.PhaseLockingData','spatialModulation','placeFields','behavior.cellinfo','ACGPeak',...
            'speedCorr','uLEDResponse.cellinfo', 'uLEDcoactivation'};
    [projectResults, projectSessionResults] = loadProjectResults('list_of_sessions',list_of_sessions,...
        'analysis_project_path', [onedrive_path 'NeuralComputationLab\ActiveProjects\ReactInh\dataCoactivation'] ,'loadLast',false, 'save_as', 'reactInh_figure_coactivation', 'list_of_results', list_of_results);
    
    % stacking pairs
    projectResults.uLEDcoactivation_pairs = stackSessionResult(projectSessionResults.uLEDcoactivation, projectSessionResults.numcells .* projectSessionResults.numcells);   
% 3. 
