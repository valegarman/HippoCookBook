%% BatchScript_analysis_SocialProject

clear; close all
targetProject= 'SocialProject';
cd('D:\FLR');
database_path = 'D:\FLR';
HCB_directory = what('SocialProject'); 

% sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
sessionsTable = readtable(['C:\Users\Jorge\Documents\GitHub\SocialProject',filesep,'indexedSessions_SocialProject.csv']);
forceReload = false;

win_resp = [-0.025 0.025];

ripple_passband = [120 200];
SW_passband = [2 10];
theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

plt = true;
force = true;

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
        ts_interSleep = [];
        count = 1;
        
        for jj = 1:length(session.epochs)
           session.epochs{jj}.behavioralParadigm
           if strcmpi(session.epochs{jj}.behavioralParadigm,'InterSleep') | strcmpi(session.epochs{jj}.behavioralParadigm,'InterSleep1') strcmpi(session.epochs{jj}.behavioralParadigm,'InterSleep2')
               ts_interSleep = [ts_interSleep; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
           end
        end
        
        if ~isempty(ts_interSleep)
            ts_interSleep = [ts_interSleep(1) ts_interSleep(end)];
        end
        

    end       
end