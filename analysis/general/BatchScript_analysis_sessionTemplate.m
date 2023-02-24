%% BatchScript_analysis_ripplesPropertiesPerSubSession
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'GLUN3Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_GLUN3Project.csv']); % the variable is called allSessions
forceReload = true;

win_resp = [-0.025 0.025];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        session = loadSession(pwd);
        
        ripples = rippleMasterDetector();
        
        
    end
        
    close all;
end