%% BatchScript_analysis_BaselinevsDrug

clear; close all
targetProject= 'MK801Project';
cd('J:\data');
database_path = 'J:\data';
HCB_directory = what('GLUN3Project'); 

% sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
sessionsTable = readtable(['C:\Users\Jorge\Documents\GitHub\MK801Project',filesep,'indexedSessions_GLUN3Project.csv']);
forceReload = false;

win_resp = [-0.025 0.025];

ripple_passband = [120 200];
SW_passband = [2 10];
theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

notchFilter = false;

plt = true;

force = true;
forceLfp = true;
forceACGPeak = false;
forceBehavior = false;
forceACG = false;

for ii = 38:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        mkdir('BaselinevsDrug')
        close all;
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
%         getAverageCCG('force',true,'skipStimulationPeriods',false);
%         spatialModulation = computeSpatialClassificationGLUN3('tint',false);
        
        ts_LongSleepPre = [];
        ts_LongSleepBaseline = [];
        ts_LongSleepDrug = [];

        count = 1;
        for jj = 1:length(session.epochs)
            session.epochs{jj}.behavioralParadigm
            if strcmpi(session.epochs{jj}.behavioralParadigm , 'PreSleep') 
                
                ts_LongSleepPre = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm , 'LongSleepBaseline') 
    
                ts_LongSleepBaseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'LongSleepDrug')
                
                ts_LongSleepDrug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
                
        % =================================
        % 2. RIPPLES 
        % =================================
        
        if isempty(dir('*.ripples_LongSleepPre.events.mat')) | isempty(dir('*.ripplesLongSleepBaseline.events.mat')) | isempty(dir('*.ripplesLongSleepDrug.events.mat')) | force
            
            % Long Sleep Pre
            if ~isempty(ts_LongSleepPre)
                try
                    targetFile = dir('*ripples.events.mat'); load(targetFile.name);
                catch
                    error('Not possible to compute ripples Baseline vs Drug.')
                end

                ts_ripples_LongSleepPre = find(InIntervals(ripples.peaks,ts_LongSleepPre));
                figure('position',[200 115 1300 800])
                plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_LongSleepPre,'inAxis',true);
                mkdir('BaselinevsDrug');
                saveas(gca,['BaselinevsDrug\plotRippleChannel_LongSleepPre.png']);


                ripples_LongSleepPre.timestamps = ripples.timestamps(ts_ripples_LongSleepPre,:);

                ripples_LongSleepPre.peaks = ripples.peaks(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.peakNormedPower = ripples.peakNormedPower(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.stdev = ripples.stdev;

                ripples_LongSleepPre.noise.times = ripples.noise.times; 
                ripples_LongSleepPre.noise.peaks = ripples.noise.peaks;
                ripples_LongSleepPre.noise.peakNormedPower = ripples.noise.peakNormedPower;

                ripples_LongSleepPre.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
                ripples_LongSleepPre.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
                ripples_LongSleepPre.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
                ripples_LongSleepPre.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
                ripples_LongSleepPre.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
                ripples_LongSleepPre.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

                ripples_LongSleepPre.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
                ripples_LongSleepPre.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;

                if isfield(ripples,'eventSpikingParameter')
                    ripples_LongSleepPre.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                    ripples_LongSleepPre.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                    ripples_LongSleepPre.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
                end

                    % Ripples stats

                ripples_LongSleepPre.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_LongSleepPre)';
                ripples_LongSleepPre.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_LongSleepPre)';
                ripples_LongSleepPre.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_LongSleepPre,:);

                ripples_LongSleepPre.rippleStats.data.incidence = length(ripples_LongSleepPre.peaks) / (ts_LongSleepPre(2)-ts_LongSleepPre(1));

                ripples_LongSleepPre.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_LongSleepPre,:),
                ripples_LongSleepPre.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_LongSleepPre);
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_LongSleepPre);
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_LongSleepPre);
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_LongSleepPre);
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_LongSleepPre,:);
                ripples_LongSleepPre.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

                corrBinSize = 0.01;

                [data,t] = CCG(ripples_LongSleepPre.peaks,ones(length(ripples_LongSleepPre.peaks),1),'binSize',corrBinSize);                    
                ripples_LongSleepPre.rippleStats.stats.acg.data = data;
                ripples_LongSleepPre.rippleStats.stats.acg.t = t;

                [rho,p] = corrcoef(ripples_LongSleepPre.rippleStats.data.peakAmplitude,ripples_LongSleepPre.rippleStats.data.peakFrequency);
                ripples_LongSleepPre.rippleStats.stats.amplitudeFrequency.rho = rho;
                ripples_LongSleepPre.rippleStats.stats.amplitudeFrequency.p = p;

                [rho,p] = corrcoef(ripples_LongSleepPre.rippleStats.data.duration,ripples_LongSleepPre.rippleStats.data.peakFrequency);
                ripples_LongSleepPre.rippleStats.stats.durationFrequency.rho = rho;
                ripples_LongSleepPre.rippleStats.stats.durationFrequency.p = p;

                [rho,p] = corrcoef(ripples_LongSleepPre.rippleStats.data.duration,ripples_LongSleepPre.rippleStats.data.peakAmplitude);
                ripples_LongSleepPre.rippleStats.stats.durationAmplitude.rho = rho;
                ripples_LongSleepPre.rippleStats.stats.durationAmplitude.p = p;

                ripples = ripples_LongSleepPre;
                save([session.general.name,'.ripples_LongSleepPre.events.mat'],'ripples');
            else
                ripples = [];
                save([session.general.name,'.ripples_LongSleepPre.events.mat'],'ripples');
            end


            % Long Sleep Baseline
            if ~isempty(ts_LongSleepBaseline)
                try
                    targetFile = dir('*ripples.events.mat'); load(targetFile.name);
                catch
                    error('Not possible to compute ripples Baseline vs Drug.')
                end

                ts_ripples_LongSleepBaseline = find(InIntervals(ripples.peaks,ts_LongSleepBaseline));

                figure('position',[200 115 1300 800])
                plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_LongSleepBaseline,'inAxis',true);
                mkdir('BaselinevsDrug');
                saveas(gca,['BaselinevsDrug\plotRippleChannel_LongSleepBaseline.png']);

                ripples_LongSleepBaseline.timestamps = ripples.timestamps(ts_ripples_LongSleepBaseline,:);

                ripples_LongSleepBaseline.peaks = ripples.peaks(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.peakNormedPower = ripples.peakNormedPower(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.stdev = ripples.stdev;

                ripples_LongSleepBaseline.noise.times = ripples.noise.times; 
                ripples_LongSleepBaseline.noise.peaks = ripples.noise.peaks;
                ripples_LongSleepBaseline.noise.peakNormedPower = ripples.noise.peakNormedPower;

                ripples_LongSleepBaseline.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
                ripples_LongSleepBaseline.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
                ripples_LongSleepBaseline.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
                ripples_LongSleepBaseline.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
                ripples_LongSleepBaseline.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
                ripples_LongSleepBaseline.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

                ripples_LongSleepBaseline.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
                ripples_LongSleepBaseline.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;

                if isfield(ripples,'eventSpikingParameters')
                    ripples_LongSleepBaseline.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                    ripples_LongSleepBaseline.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                    ripples_LongSleepBaseline.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
                end

                    % Ripples stats

                ripples_LongSleepBaseline.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_LongSleepBaseline)';
                ripples_LongSleepBaseline.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_LongSleepBaseline)';
                ripples_LongSleepBaseline.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.data.incidence = length(ripples_LongSleepBaseline.peaks) / (ts_LongSleepBaseline(2)-ts_LongSleepBaseline(1));

                ripples_LongSleepBaseline.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_LongSleepBaseline,:),
                ripples_LongSleepBaseline.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_LongSleepBaseline);
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_LongSleepBaseline);
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_LongSleepBaseline);
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_LongSleepBaseline);
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_LongSleepBaseline,:);
                ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

                corrBinSize = 0.01;

                [data,t] = CCG(ripples_LongSleepBaseline.peaks,ones(length(ripples_LongSleepBaseline.peaks),1),'binSize',corrBinSize);                    
                ripples_LongSleepBaseline.rippleStats.stats.acg.data = data;
                ripples_LongSleepBaseline.rippleStats.stats.acg.t = t;

                [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.peakAmplitude,ripples_LongSleepBaseline.rippleStats.data.peakFrequency);
                ripples_LongSleepBaseline.rippleStats.stats.amplitudeFrequency.rho = rho;
                ripples_LongSleepBaseline.rippleStats.stats.amplitudeFrequency.p = p;

                [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.duration,ripples_LongSleepBaseline.rippleStats.data.peakFrequency);
                ripples_LongSleepBaseline.rippleStats.stats.durationFrequency.rho = rho;
                ripples_LongSleepBaseline.rippleStats.stats.durationFrequency.p = p;

                [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.duration,ripples_LongSleepBaseline.rippleStats.data.peakAmplitude);
                ripples_LongSleepBaseline.rippleStats.stats.durationAmplitude.rho = rho;
                ripples_LongSleepBaseline.rippleStats.stats.durationAmplitude.p = p;

                ripples = ripples_LongSleepBaseline;
                save([session.general.name,'.ripples_LongSleepBaseline.events.mat'],'ripples');
            else
                ripples = [];
                save([session.general.name,'.ripples_LongSleepBaseline.events.mat'],'ripples');
            end
        end
               
        % Long Sleep Drug
        if ~isempty(ts_LongSleepDrug)
            try
                targetFile = dir('*ripples.events.mat'); load(targetFile.name);
            catch
                error('Not possible to compute ripples Baseline vs Drug.')
            end

            ts_ripples_LongSleepDrug = find(InIntervals(ripples.peaks,ts_LongSleepDrug));

            figure('position',[200 115 1300 800])
            plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_LongSleepDrug,'inAxis',true);
            mkdir('BaselinevsDrug');
            saveas(gca,['BaselinevsDrug\plotRippleChannel_LongSleepDrug.png']);

            ripples_LongSleepDrug.timestamps = ripples.timestamps(ts_ripples_LongSleepDrug,:);

            ripples_LongSleepDrug.peaks = ripples.peaks(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.peakNormedPower = ripples.peakNormedPower(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.stdev = ripples.stdev;

            ripples_LongSleepDrug.noise.times = ripples.noise.times; 
            ripples_LongSleepDrug.noise.peaks = ripples.noise.peaks;
            ripples_LongSleepDrug.noise.peakNormedPower = ripples.noise.peakNormedPower;

            ripples_LongSleepDrug.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
            ripples_LongSleepDrug.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
            ripples_LongSleepDrug.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
            ripples_LongSleepDrug.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
            ripples_LongSleepDrug.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
            ripples_LongSleepDrug.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

            ripples_LongSleepDrug.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
            ripples_LongSleepDrug.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;

            if isfield(ripples,'eventSpikingParameters')
                ripples_LongSleepDrug.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                ripples_LongSleepDrug.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                ripples_LongSleepDrug.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
            end

                % Ripples stats

            ripples_LongSleepDrug.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_LongSleepDrug)';
            ripples_LongSleepDrug.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_LongSleepDrug)';
            ripples_LongSleepDrug.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.data.incidence = length(ripples_LongSleepDrug.peaks) / (ts_LongSleepDrug(2)-ts_LongSleepDrug(1));

            ripples_LongSleepDrug.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_LongSleepDrug,:),
            ripples_LongSleepDrug.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_LongSleepDrug);
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_LongSleepDrug);
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_LongSleepDrug);
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_LongSleepDrug);
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_LongSleepDrug,:);
            ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

            corrBinSize = 0.01;

            [data,t] = CCG(ripples_LongSleepDrug.peaks,ones(length(ripples_LongSleepDrug.peaks),1),'binSize',corrBinSize);                    
            ripples_LongSleepDrug.rippleStats.stats.acg.data = data;
            ripples_LongSleepDrug.rippleStats.stats.acg.t = t;

            [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.peakAmplitude,ripples_LongSleepDrug.rippleStats.data.peakFrequency);
            ripples_LongSleepDrug.rippleStats.stats.amplitudeFrequency.rho = rho;
            ripples_LongSleepDrug.rippleStats.stats.amplitudeFrequency.p = p;

            [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.duration,ripples_LongSleepDrug.rippleStats.data.peakFrequency);
            ripples_LongSleepDrug.rippleStats.stats.durationFrequency.rho = rho;
            ripples_LongSleepDrug.rippleStats.stats.durationFrequency.p = p;

            [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.duration,ripples_LongSleepDrug.rippleStats.data.peakAmplitude);
            ripples_LongSleepDrug.rippleStats.stats.durationAmplitude.rho = rho;
            ripples_LongSleepDrug.rippleStats.stats.durationAmplitude.p = p;

            ripples = ripples_LongSleepDrug;
            save([session.general.name,'.ripples_LongSleepDrug.events.mat'],'ripples');
        else
            ripples = [];
            save([session.general.name,'.ripples_LongSleepDrug.events.mat'],'ripples');
        end
    end
        
        % =================================
        % 3. RIPPLES psth 
        % =================================
        
        if isempty(dir('*.ripples_LongSleepPre_psth.cellinfo.mat')) | isempty(dir('*.ripples_LongSleepBaseline_psth.cellinfo.mat')) | isempty(dir('*.ripples_LongSleepDrug_psth.cellinfo.mat')) | force
            
            ripples = rippleMasterDetector('rippleChannel',[],'SWChannel',[],'force',false,'removeOptogeneticStimulation',true);
            
            % Long Sleep Pre
    
            psthRipples_LongSleepPre = spikesPsth([],'eventType','ripples','restrictIntervals',ts_LongSleepPre,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

            t = psthRipples_LongSleepPre.timestamps;
            winSizePlot = [-0.5 0.5];
            figure('position',[200 115 1300 800])
            subplot(1,2,1)
            imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepPre.responsecurve,1)],...
                psthRipples_LongSleepPre.responsecurveSmooth); caxis([0 10]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');

            subplot(1,2,2)
            imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepPre.responsecurve,1)],...
                psthRipples_LongSleepPre.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');
            mkdir('BaselineVsDrug');

            saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_LongSleepPre.png']);

            raster = psthRipples_LongSleepPre.raster;
            psth = rmfield(psthRipples_LongSleepPre,'raster');

            save([session.general.name,'.ripples_LongSleepPre_psth.cellinfo.mat'],'psth','-v7.3');
            save([session.general.name,'.ripples_LongSleepPre_raster.cellinfo.mat'],'raster','-v7.3');

            % Long Sleep Baseline
            if ~isempty(ts_LongSleepBaseline)
                psthRipples_LongSleepBaseline = spikesPsth([],'eventType','ripples','restrictIntervals',ts_LongSleepBaseline,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

                t = psthRipples_LongSleepBaseline.timestamps;
                winSizePlot = [-0.5 0.5];
                figure('position',[200 115 1300 800])
                subplot(1,2,1)
                imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepBaseline.responsecurve,1)],...
                    psthRipples_LongSleepBaseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');

                subplot(1,2,2)
                imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepBaseline.responsecurve,1)],...
                    psthRipples_LongSleepBaseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');
                mkdir('BaselineVsDrug');

                saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_LongSleepBaseline.png']);

                raster = psthRipples_LongSleepBaseline.raster;
                psth = rmfield(psthRipples_LongSleepBaseline,'raster');

                save([session.general.name,'.ripples_LongSleepBaseline_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_LongSleepBaseline_raster.cellinfo.mat'],'raster','-v7.3');
            else
                raster = [];
                psth = [];

                save([session.general.name,'.ripples_LongSleepBaseline_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_LongSleepBaseline_raster.cellinfo.mat'],'raster','-v7.3');
            end
            
            
            % Long Sleep Drug
            if ~isempty(ts_LongSleepDrug)
                psthRipples_LongSleepDrug = spikesPsth([],'eventType','ripples','restrictIntervals',ts_LongSleepDrug,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

                t = psthRipples_LongSleepDrug.timestamps;
                winSizePlot = [-0.5 0.5];
                figure('position',[200 115 1300 800])
                subplot(1,2,1)
                imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepDrug.responsecurve,1)],...
                    psthRipples_LongSleepDrug.responsecurveSmooth); caxis([0 10]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');

                subplot(1,2,2)
                imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepDrug.responsecurve,1)],...
                    psthRipples_LongSleepDrug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');
                mkdir('BaselineVsDrug');

                saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_LongSleepBaseline.png']);

                raster = psthRipples_LongSleepDrug.raster;
                psth = rmfield(psthRipples_LongSleepDrug,'raster');

                save([session.general.name,'.ripples_LongSleepDrug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_LongSleepDrug_raster.cellinfo.mat'],'raster','-v7.3');
            else
                raster = [];
                psth = [];

                save([session.general.name,'.ripples_LongSleepDrug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_LongSleepDrug_raster.cellinfo.mat'],'raster','-v7.3');
            end
        end
        
        
    close all;
    clc;
end