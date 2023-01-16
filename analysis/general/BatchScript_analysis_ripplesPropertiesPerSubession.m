%% BatchScript_analysis_ripplesPropertiesPerSubSession
% place your code to run an analysis across all sessions for a given
% project

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
        
        if isempty(dir('*.ripplesSubsessions.events.mat')) || forceReload
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
                    
                    ripplesSubsessions.(session.epochs{jj}.name).timestamps = ripples.timestamps(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).peaks = ripples.peaks(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).peakNormedPower = ripples.peakNormedPower(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).stdev = ripples.stdev;
                    
                    ripplesSubsessions.(session.epochs{jj}.name).noise.times = ripples.noise.times; 
                    ripplesSubsessions.(session.epochs{jj}.name).noise.peaks = ripples.noise.peaks;
                    ripplesSubsessions.(session.epochs{jj}.name).noise.peakNormedPower = ripples.noise.peakNormedPower;
                    
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.detectorname = ripples.detectorinfo.detectorname;
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
                    ripplesSubsessions.(session.epochs{jj}.name).detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
                    
                    ripplesSubsessions.(session.epochs{jj}.name).artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
                    ripplesSubsessions.(session.epochs{jj}.name).artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
                    
                    ripplesSubsessions.(session.epochs{jj}.name).eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                    ripplesSubsessions.(session.epochs{jj}.name).eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                    ripplesSubsessions.(session.epochs{jj}.name).eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
                    
                    % Ripple Stats
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples)';
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples)';
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples,:);
                    
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples,:),
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
                    
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples,:);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
                    
                    % This needs to be computed for each subsession
                    corrBinSize = 0.01;
                    
                    [data,t] = CCG(ripplesSubsessions.(session.epochs{jj}.name).peaks,ones(length(ripplesSubsessions.(session.epochs{jj}.name).peaks),1),'binSize',corrBinSize);                    
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.acg.data = data;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.acg.t = t;
                    
                    [rho,p] = corrcoef(ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakAmplitude,ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakFrequency);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.amplitudeFrequency.rho = rho;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.amplitudeFrequency.p = p;
                    
                    [rho,p] = corrcoef(ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.duration,ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakFrequency);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.durationFrequency.rho = rho;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.durationFrequency.p = p;
                    
                    [rho,p] = corrcoef(ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.duration,ripplesSubsessions.(session.epochs{jj}.name).rippleStats.data.peakAmplitude);
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.durationAmplitude.rho = rho;
                    ripplesSubsessions.(session.epochs{jj}.name).rippleStats.stats.durationAmplitude.p = p;
                    
                    
                catch
                    warning('Not enough ripples in this epoch');
                    ripplesSubsessions.(session.epochs{jj}.name) = [];
                end
            end

            save([session.general.name, '.ripplesSubsessions.events.mat'],'ripplesSubsessions','-v7.3');
            
            clear ripplesSubsessions

            catch
                warning('Analysis was not possible!');
            end
        end
    end
    close all;
end
%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions