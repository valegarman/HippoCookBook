%% BatchScript_analysis_ripplesPerSubsession

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'POTometry';

for ii = 95:length(sessionsTable.SessionName)
    if contains(sessionsTable.Project(ii), targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}]);

        % 
        clear basepath session basename cell_metrics spikes
        basepath = pwd;
        session = loadSession;
        basename = basenameFromBasepath(basepath);
        cell_metrics = loadCellMetrics;
        parameters.forceReload = true;
        spikes = loadSpikes;
        
        % Compute unit classification

        
    end
    close all;
end