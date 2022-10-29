%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'All';

HCB_directory = what('HippoCookBook'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        try

            %%% your code goes here...
            getSpikesRank('events','upstates')
            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end

%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions