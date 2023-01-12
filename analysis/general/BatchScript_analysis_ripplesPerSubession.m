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
        
        if isempty(dir('*.ripples_psthSubsessions.cellinfo.mat')) || forceReload
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
%                     plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples); 
                    psthRipples = spikesPsth(ripples.peaks(ts_ripples),'eventType','ripples','min_pulsesNumber',0,'force',true,'saveMat',false,'savePlot',false);
                    
                    ripples_timestamps = psthRipples.timestamps;
                    win_Z = find(ripples_timestamps<=-0.1);
                    
                    % Computing responseZ (instead of doing in master
                    % function due to subsession analysis)
                    for kk = 1:size(psthRipples.responsecurveSmooth,1)
                        psthRipples.responseZ(kk,:) = (psthRipples.responsecurveSmooth(kk,:) - ...
                            mean(psthRipples.responsecurveSmooth(kk,win_Z)))./std(psthRipples.responsecurveSmooth(kk,win_Z));
                    end
                    
                    % Computing peakResponse and peakResponseZ (instead of
                    % doing in master function due to subsession analysis)
                    win = find(ripples_timestamps>=win_resp(1) & ripples_timestamps<=win_resp(2));
                    for kk = 1:size(psthRipples.responsecurve,1)
                        psthRipples.peakResponse(kk,:) = nanmean(psthRipples.responsecurve(kk,win),2); % delta peak response
                    end
                    for kk = 1:size(psthRipples.responsecurveZ,1)
                        psthRipples.peakResponseZ(kk,:) = nanmean(psthRipples.responseZ(kk,win),2); % delta peak response
                    end
                    % Saving in subfolder
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
            
            clear psthRipplesSubsessions
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