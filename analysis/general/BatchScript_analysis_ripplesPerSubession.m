%% BatchScript_analysis_ripplesPerSubsession
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        if isempty(dir('*.ripples_psthSubsessions.cellinfo.mat'))
            try
            %%% your code goes here...
            % Load ripples
            targetFile = dir('*.ripples.events.mat');
            try
                load(targetFile.name);
            catch
                warning('ripples.events.mat not found. Quitting analysis...');
            end
            % Load session info
            session = loadSession();

            for jj = 1:length(session.epochs)
                ts = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];                
                ts_ripples = find(InIntervals(ripples.peaks,ts));
                try
                    plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples); 
                    psthRipples = spikesPsth(ripples.peaks(ts_ripples),'eventType','ripples','min_pulsesNumber',0,'force',true,'saveMat',false,'savePlot',false);
                    psthRipplesSubsessions.(session.epochs{jj}.name) = psthRipples;
                    psthrasterRipplesSubsessions.(session.epochs{jj}.name) = psthRipples.raster;
                catch
                    warning('Not enough ripples in this epoch');
                    psthRipples = [];
                    psthRipples.raster = [];
                    psthRipplesSubsessions.(session.epochs{jj}.name) = psthRipples;
                    psthrasterRipplesSubsessions.(session.epochs{jj}.name) = psthRipples.raster;
                end
            end

            save([session.general.name, '.ripples_psthSubsessions.cellinfo.mat'],'psthRipplesSubsessions','-v7.3');
            save([session.general.name, '.ripples_rasterSubsessions.cellinfo.mat'],'psthrasterRipplesSubsessions','-v7.3'); 
            
            clear psthRipplesSubsesions
            clear psthrasterRipplesSubsessions

            catch
                warning('Analysis was not possible!');
            end
        end
    end
    close all;
end
%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions