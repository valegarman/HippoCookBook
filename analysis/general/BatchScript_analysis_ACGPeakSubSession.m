%% BatchScript_analysis_ripplesPerSubsession

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
forceReload = false;

win_resp = [-0.025 0.025];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        if isempty(dir('*.ACGPeak_SubSession.cellinfo.mat')) | forceReload
            try
                % ACGPeak SubSessions
                acgPeak = getACGPeakSubSession;
            catch
                warning('Not possible to run getACGPeakSubSession...');
            end
            
        end
    end
    close all;
end