%% BatchScript_analysis_ripplesPropertiesPerSubSession
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
forceReload = true;

win_resp = [-0.025 0.025];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
%         session = sessionTemplate(pwd,'showGUI',true); 
        session = loadSession(pwd);
        if length(session.epochs) == 10
            for jj = 1:length(session.epochs)
                switch jj
                    case 1
                        session.epochs{jj}.behavioralParadigm = 'PreSleep';
                    case 2
                        session.epochs{jj}.behavioralParadigm = 'Maze1Baseline';
                    case 3
                        session.epochs{jj}.behavioralParadigm = 'InterMazeBaseline';
                    case 4
                        session.epochs{jj}.behavioralParadigm = 'Maze2Baseline';
                    case 5
                        session.epochs{jj}.behavioralParadigm = 'LongSleepBaseline';
                    case 6
                        session.epochs{jj}.behavioralParadigm = 'InjectionInterTrial';
                    case 7
                        session.epochs{jj}.behavioralParadigm = 'Maze1Drug';
                    case 8
                        session.epochs{jj}.behavioralParadigm = 'InterMazeDrug';
                    case 9
                        session.epochs{jj}.behavioralParadigm = 'Maze2Drug';
                    case 10
                        session.epochs{jj}.behavioralParadigm = 'LongSleepDrug';
                end
            end
        elseif length(session.epochs) == 5 & ~strcmpi(session.general.name,'IPO429_131021_sess2')
            for jj = 1:length(session.epochs)
                switch jj
                    case 1
                        session.epochs{jj}.behavioralParadigm = 'PreSleep';
                    case 2
                        session.epochs{jj}.behavioralParadigm = 'Maze1Baseline';
                    case 3
                        session.epochs{jj}.behavioralParadigm = 'InterMazeBaseline';
                    case 4
                        session.epochs{jj}.behavioralParadigm = 'Maze2Baseline';
                    case 5
                        session.epochs{jj}.behavioralParadigm = 'LongSleepBaseline';
                end
            end
        elseif length(session.epochs) == 5 & strcmpi(session.general.name,'IPO429_131021_sess2')
            for jj = 1:length(session.epochs)
                switch jj
                    case 1
                        session.epochs{jj}.behavioralParadigm = 'PreSleep';
                    case 2
                        session.epochs{jj}.behavioralParadigm = 'Maze1Baseline';
                    case 3
                        session.epochs{jj}.behavioralParadigm = 'LongSleepBaseline';
                    case 4
                        session.epochs{jj}.behavioralParadigm = 'Maze1Drug';
                    case 5
                        session.epochs{jj}.behavioralParadigm = 'LongSleepDrug';
                end
            end
        elseif length(session.epochs) == 6 
            for jj = 1:length(session.epochs)
                switch jj
                    case 1
                        session.epochs{jj}.behavioralParadigm = 'PreSleep';
                    case 2
                        session.epochs{jj}.behavioralParadigm = 'Maze1Baseline';
                    case 3
                        session.epochs{jj}.behavioralParadigm = 'LongSleepBaseline';
                    case 4
                        session.epochs{jj}.behavioralParadigm = 'InjectionInterTrial';
                    case 5
                        session.epochs{jj}.behavioralParadigm = 'Maze1Drug';
                    case 6
                        session.epochs{jj}.behavioralParadigm = 'LongSleepDrug';
                end
            end
        end
    end
    save([session.general.name,'.session.mat'],'session');
    close all;
end