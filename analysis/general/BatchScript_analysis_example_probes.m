
%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'All';

HCB_directory = what('HippoCookBook'); 

load([HCB_directory.path filesep 'indexedSessions.mat']);
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        try

            %%% your code goes here...
            session = loadSession;
            try
                if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
                    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
                    load(file.name);
                end
                    excludeManipulationIntervals = optoPulses.stimulationEpochs;
            catch
                warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
            end
            cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals);
            close all
            clear session
            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end