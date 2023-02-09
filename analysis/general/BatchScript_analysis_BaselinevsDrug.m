%% BatchScript_analysis_BaselinevsDrug

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
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

plt = true;

force = true;

for ii = 5:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        ts_Baseline = [];
        ts_Drug = [];
        count = 1;
        for jj = 1:length(session.epochs)
            session.epochs{jj}.behavioralParadigm
            if strcmpi(session.epochs{jj}.behavioralParadigm , 'PreSleep') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze1Baseline')  | strcmpi(session.epochs{jj}.behavioralParadigm, 'InterMazeBaseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze2Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze3Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'LongSleepBaseline')
                
                ts_Baseline = [ts_Baseline ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
%             elseif strcmpi(session.epochs{jj}.behavioralParadigm , 'InjectionInterTrial') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze1Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'InterMazeDrug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze2Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze3Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'LongSleepDrug')
            elseif strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze1Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'InterMazeDrug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze2Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze3Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'LongSleepDrug')
    
                ts_Drug = [ts_Drug session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        ts_Baseline = [ts_Baseline(1) ts_Baseline(end)];
        if ~isempty(ts_Drug)
            ts_Drug = [ts_Drug(1) ts_Drug(end)];
        end
        
        % ======================================
        % 0. SESSION PERFORMANCE
        % ======================================
        
        % Baseline
        try
             if isempty(dir('*.OpenField_Baseline.mat')) | isempty(dir('*.OpenField_Drug.mat'))
                performance = getSessionPerformance('includeIntervals',ts_Baseline);

                OpenField = performance.OpenField;
                if isfield(performance,'YMaze')
                    YMaze = performance.YMaze;
                end
                performance = OpenField;
                save([session.general.name,'.OpenField_Baseline.mat'],'performance');

                if exist('YMaze','var')  
                    performance = YMaze;
                    save([session.general.name,'.YMaze_Baseline.mat'],'performance');
                end
            end

            % Drug
           if isempty(dir('*.YMaze_Baseline.mat')) | isempty(dir('*.YMaze_Drug.mat')) 
               performance = getSessionPerformance('includeIntervals',ts_Drug);

                OpenField = performance.OpenField;
                if isfield(performance,'YMaze')
                    YMaze = performance.YMaze;
                end
                performance = OpenField;
                save([session.general.name,'.OpenField_Drug.mat'],'performance');

                if exist('YMaze','var')  
                    performance = YMaze;
                    save([session.general.name,'.YMaze_Drug.mat'],'performance');
                end
           end
        catch
        end
        
        
        % ======================================
        % 1. UP-DOWN psth
        % ======================================
        
        if isempty(dir('*.slowOscillations_Baseline_psth.cellinfo.mat')) | isempty(dir('*.slowOscillations_Drug_psth.cellinfo.mat')) | force
            
            UDStates = detectUD('plotOpt', true,'forceDetect',false','NREMInts','all');

            % Baseline

            psthUD_Baseline = spikesPsth([],'eventType','slowOscillations','restrictIntervals',ts_Baseline,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

            t = psthUD_Baseline.timestamps;
            winSizePlot = [-0.5 0.5];
            figure,
            subplot(1,2,1)
            imagesc([t(1) t(end)],[1 size(psthUD_Baseline.responsecurve,1)],...
                psthUD_Baseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');

            subplot(1,2,2)
            imagesc([t(1) t(end)],[1 size(psthUD_Baseline.responsecurve,1)],...
                psthUD_Baseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');
            mkdir('BaselineVsDrug');

            saveas(gcf,['BaselinevsDrug\spikesPsthRate_slowOscillations_Baseline.png']);

            raster = psthUD_Baseline.raster;
            psth = rmfield(psthUD_Baseline,'raster');

            save([session.general.name,'.slowOscillations_Baseline_psth.cellinfo.mat'],'psth','-v7.3');
            save([session.general.name,'.slowOscillations_Baseline_raster.cellinfo.mat'],'raster','-v7.3');
            
            
            % Drug
            if ~isempty(ts_Drug)
                psthUD_Drug = spikesPsth([],'eventType','slowOscillations','restrictIntervals',ts_Drug,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

                t = psthUD_Drug.timestamps;
                winSizePlot = [-0.5 0.5];
                figure,
                subplot(1,2,1)
                imagesc([t(1) t(end)],[1 size(psthUD_Drug.responsecurve,1)],...
                    psthUD_Drug.responsecurveSmooth); caxis([0 10]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');

                subplot(1,2,2)
                imagesc([t(1) t(end)],[1 size(psthUD_Drug.responsecurve,1)],...
                    psthUD_Drug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');
                mkdir('BaselineVsDrug');

                saveas(gcf,['BaselinevsDrug\spikesPsthRate_slowOscillations_Druge.png']);

                raster = psthUD_Drug.raster;
                psth = rmfield(psthUD_Drug,'raster');

                save([session.general.name,'.slowOscillations_Drug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.slowOscillations_Drug_raster.cellinfo.mat'],'raster','-v7.3');
            else
                raster = [];
                psth = [];
                save([session.general.name,'.slowOscillations_Drug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.slowOscillations_Drug_raster.cellinfo.mat'],'raster','-v7.3');
            end
        end
        
        % =================================
        % 2. RIPPLES 
        % =================================
        
        if isempty(dir('*.ripples_Baseline.events.mat')) | isempty(dir('*.ripples_Drug.events.mat')) | force
            
            % Baseline
            try
                targetFile = dir('*ripples.events.mat'); load(targetFile.name);
            catch
                error('Not possible to compute ripples Baseline vs Drug.')
            end

            ts_ripples_Baseline = find(InIntervals(ripples.peaks,ts_Baseline));
            figure('position',[200 115 1300 800])
            plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_Baseline,'inAxis',true);
            mkdir('BaselinevsDrug');
            saveas(gcf,['BaselinevsDrug\plotRippleChannel_Baseline.png']);


            ripples_Baseline.timestamps = ripples.timestamps(ts_ripples_Baseline,:);

            ripples_Baseline.peaks = ripples.peaks(ts_ripples_Baseline,:);
            ripples_Baseline.peakNormedPower = ripples.peakNormedPower(ts_ripples_Baseline,:);
            ripples_Baseline.stdev = ripples.stdev;

            ripples_Baseline.noise.times = ripples.noise.times; 
            ripples_Baseline.noise.peaks = ripples.noise.peaks;
            ripples_Baseline.noise.peakNormedPower = ripples.noise.peakNormedPower;

            ripples_Baseline.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
            ripples_Baseline.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
            ripples_Baseline.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
            ripples_Baseline.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
            ripples_Baseline.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
            ripples_Baseline.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

            ripples_Baseline.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
            ripples_Baseline.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
            
            if isfield(ripples,'eventSpikingParameter')
                ripples_Baseline.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                ripples_Baseline.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                ripples_Baseline.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
            end

                % Ripples stats

            ripples_Baseline.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_Baseline)';
            ripples_Baseline.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_Baseline)';
            ripples_Baseline.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_Baseline,:);
            
            ripples_Baseline.rippleStats.data.incidence = length(ripples_Baseline.peaks) / (ts_Baseline(2)-ts_Baseline(1));

            ripples_Baseline.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_Baseline,:),
            ripples_Baseline.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

            ripples_Baseline.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_Baseline);
            ripples_Baseline.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
            ripples_Baseline.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
            ripples_Baseline.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
            ripples_Baseline.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_Baseline);
            ripples_Baseline.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_Baseline);
            ripples_Baseline.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_Baseline);
            ripples_Baseline.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_Baseline,:);
            ripples_Baseline.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

            corrBinSize = 0.01;

            [data,t] = CCG(ripples_Baseline.peaks,ones(length(ripples_Baseline.peaks),1),'binSize',corrBinSize);                    
            ripples_Baseline.rippleStats.stats.acg.data = data;
            ripples_Baseline.rippleStats.stats.acg.t = t;

            [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.peakAmplitude,ripples_Baseline.rippleStats.data.peakFrequency);
            ripples_Baseline.rippleStats.stats.amplitudeFrequency.rho = rho;
            ripples_Baseline.rippleStats.stats.amplitudeFrequency.p = p;

            [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.duration,ripples_Baseline.rippleStats.data.peakFrequency);
            ripples_Baseline.rippleStats.stats.durationFrequency.rho = rho;
            ripples_Baseline.rippleStats.stats.durationFrequency.p = p;

            [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.duration,ripples_Baseline.rippleStats.data.peakAmplitude);
            ripples_Baseline.rippleStats.stats.durationAmplitude.rho = rho;
            ripples_Baseline.rippleStats.stats.durationAmplitude.p = p;

            ripples = ripples_Baseline;
            save([session.general.name,'.ripples_Baseline.events.mat'],'ripples');


            % Drug
            if ~isempty(ts_Drug)
                try
                    targetFile = dir('*ripples.events.mat'); load(targetFile.name);
                catch
                    error('Not possible to compute ripples Baseline vs Drug.')
                end

                ts_ripples_Drug = find(InIntervals(ripples.peaks,ts_Drug));

                figure('position',[200 115 1300 800])
                plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_Drug,'inAxis',true);
                mkdir('BaselinevsDrug');
                saveas(gcf,['BaselinevsDrug\plotRippleChannel_Drug.png']);

                ripples_Drug.timestamps = ripples.timestamps(ts_ripples_Drug,:);

                ripples_Drug.peaks = ripples.peaks(ts_ripples_Drug,:);
                ripples_Drug.peakNormedPower = ripples.peakNormedPower(ts_ripples_Drug,:);
                ripples_Drug.stdev = ripples.stdev;

                ripples_Drug.noise.times = ripples.noise.times; 
                ripples_Drug.noise.peaks = ripples.noise.peaks;
                ripples_Drug.noise.peakNormedPower = ripples.noise.peakNormedPower;

                ripples_Drug.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
                ripples_Drug.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
                ripples_Drug.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
                ripples_Drug.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
                ripples_Drug.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
                ripples_Drug.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

                ripples_Drug.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
                ripples_Drug.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;

                if isfield(ripples,'eventSpikingParameters')
                    ripples_Drug.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                    ripples_Drug.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                    ripples_Drug.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
                end

                    % Ripples stats

                ripples_Drug.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_Drug)';
                ripples_Drug.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_Drug)';
                ripples_Drug.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.data.incidence = length(ripples_Drug.peaks) / (ts_Drug(2)-ts_Drug(1));

                ripples_Drug.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_Drug,:),
                ripples_Drug.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

                ripples_Drug.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_Drug);
                ripples_Drug.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
                ripples_Drug.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
                ripples_Drug.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
                ripples_Drug.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_Drug);
                ripples_Drug.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_Drug);
                ripples_Drug.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_Drug);
                ripples_Drug.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_Drug,:);
                ripples_Drug.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

                corrBinSize = 0.01;

                [data,t] = CCG(ripples_Drug.peaks,ones(length(ripples_Drug.peaks),1),'binSize',corrBinSize);                    
                ripples_Drug.rippleStats.stats.acg.data = data;
                ripples_Drug.rippleStats.stats.acg.t = t;

                [rho,p] = corrcoef(ripples_Drug.rippleStats.data.peakAmplitude,ripples_Drug.rippleStats.data.peakFrequency);
                ripples_Drug.rippleStats.stats.amplitudeFrequency.rho = rho;
                ripples_Drug.rippleStats.stats.amplitudeFrequency.p = p;

                [rho,p] = corrcoef(ripples_Drug.rippleStats.data.duration,ripples_Drug.rippleStats.data.peakFrequency);
                ripples_Drug.rippleStats.stats.durationFrequency.rho = rho;
                ripples_Drug.rippleStats.stats.durationFrequency.p = p;

                [rho,p] = corrcoef(ripples_Drug.rippleStats.data.duration,ripples_Drug.rippleStats.data.peakAmplitude);
                ripples_Drug.rippleStats.stats.durationAmplitude.rho = rho;
                ripples_Drug.rippleStats.stats.durationAmplitude.p = p;

                ripples = ripples_Drug;
                save([session.general.name,'.ripples_Drug.events.mat'],'ripples');
            else
                ripples = [];
                save([session.general.name,'.ripples_Drug.events.mat'],'ripples');
            end
        end
               
        % =================================
        % 3. RIPPLES psth 
        % =================================
        
        if isempty(dir('*.ripples_Baseline_psth.cellinfo.mat')) | isempty(dir('*.ripples_Drug_psth.cellinfo.mat')) | force
            
            ripples = rippleMasterDetector('rippleChannel',[],'SWChannel',[],'force',false,'removeOptogeneticStimulation',true);
            
            % Baseline
    
            psthRipples_Baseline = spikesPsth([],'eventType','ripples','restrictIntervals',ts_Baseline,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

            t = psthRipples_Baseline.timestamps;
            winSizePlot = [-0.5 0.5];
            figure,
            subplot(1,2,1)
            imagesc([t(1) t(end)],[1 size(psthRipples_Baseline.responsecurve,1)],...
                psthRipples_Baseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');

            subplot(1,2,2)
            imagesc([t(1) t(end)],[1 size(psthRipples_Baseline.responsecurve,1)],...
                psthRipples_Baseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
            xlabel('Time');
            mkdir('BaselineVsDrug');

            saveas(gcf,['BaselinevsDrug\spikesPsthRate_ripples_Baseline.png']);

            raster = psthRipples_Baseline.raster;
            psth = rmfield(psthRipples_Baseline,'raster');

            save([session.general.name,'.ripples_Baseline_psth.cellinfo.mat'],'psth','-v7.3');
            save([session.general.name,'.ripples_Baseline_raster.cellinfo.mat'],'raster','-v7.3');

            % Drug 
            if ~isempty(ts_Drug)
                psthRipples_Drug = spikesPsth([],'eventType','ripples','restrictIntervals',ts_Drug,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

                t = psthRipples_Drug.timestamps;
                winSizePlot = [-0.5 0.5];
                figure,
                subplot(1,2,1)
                imagesc([t(1) t(end)],[1 size(psthRipples_Drug.responsecurve,1)],...
                    psthRipples_Drug.responsecurveSmooth); caxis([0 10]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');

                subplot(1,2,2)
                imagesc([t(1) t(end)],[1 size(psthRipples_Drug.responsecurve,1)],...
                    psthRipples_Drug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
                ylabel('Cells');
                xlabel('Time');
                mkdir('BaselineVsDrug');

                saveas(gcf,['BaselinevsDrug\spikesPsthRate_ripples_Drug.png']);

                raster = psthRipples_Drug.raster;
                psth = rmfield(psthRipples_Drug,'raster');

                save([session.general.name,'.ripples_Drug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_Drug_raster.cellinfo.mat'],'raster','-v7.3');
            else
                raster = [];
                psth = [];

                save([session.general.name,'.ripples_Drug_psth.cellinfo.mat'],'psth','-v7.3');
                save([session.general.name,'.ripples_Drug_raster.cellinfo.mat'],'raster','-v7.3');
            end
        end
        
        % =================================
        % 4. PHASE MODULATION
        % =================================
        
        if isempty(dir('*.theta_6-12_Baseline.PhaseLockingData.cellinfo.mat')) | isempty(dir('*.theta_6-12_Drug.PhaseLockingData.cellinfo.mat')) | force
            try
                targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
                rippleChannel = session.analysisTags.rippleChannel;
                if rippleChannel ~= ripples.detectorinfo.detectionchannel;
                    error('session metadata and ripple channels are not the same');
                end
            catch
                rippleChannel = ripples.detectorinfo.detectionchannel;
            end
            
            try
                targetFile = dir('*.sharpwaves.events.mat'); load(targetFile.name);
                if ~isempty(targetFile)
                    SWChannel = session.analysisTags.SWChannel;
                    if SWChannel ~= SW.detectorinfo.detectionchannel
                        error('session metadata and ripple channels are not the same');
                    end
                else
                    SWChannel = [];
                end
            catch
                SWChannel = [];
            end
            
            try
                targetFile = dir('*.thetaEpochs.states.mat'); load(targetFile.name);
                thetaChannel = session.analysisTags.thetaChannel;
                if thetaChannel ~= thetaEpochs.channel
                    error('session metadata and theta channels are not the same');
                end
            catch
                thetaChannel = thetaEpochs.channel;
            end
            
            
            % BASELINE
            
            [phaseMod_Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Baseline,'plotting',false,'saveMat',false);

            rippleMod = phaseMod_Baseline.ripples;
            SWMod = phaseMod_Baseline.SharpWave;
            thetaMod = phaseMod_Baseline.theta;
            lgammaMod = phaseMod_Baseline.lgamma;
            hgammaMod = phaseMod_Baseline.hgamma;
            thetaRunMod = phaseMod_Baseline.thetaRunMod;
            thetaREMMod = phaseMod_Baseline.thetaREMMod;

            save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
            
            if plt
                
                % Ripples modulation
               try
                   figure,
                   set(gcf,'Position',get(0,'ScreenSize'))
                   for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                   end 
                   saveas(gcf,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run ripples modulation Baseline vs Drug...');
               end
            end
            
            save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
            
            if plt
                % SW Modulation
               try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run SW modulation Baseline vs Drug');
               end
            end

            save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
            
            if plt
                % Theta modulation
               try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
            if plt
                % lgamma modulation

               try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Baseline_PhaseModulation.png'])
               catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
            if plt
                % hgamma modulation

               try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
            
            % ThetaRun modulation
            if plt
               try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run thetaRun modulation Baseline vs Drug...');
               end
            end
        

            save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
            
            if plt
                % ThetaREM modulation
                try
                    figure,
                    set(gcf,'Position',get(0,'ScreenSize'))
                    for i = 1:length(spikes.UID)
                        subplot(7,ceil(size(spikes.UID,2)/7),i)
                        area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
                        hold on;
                        ax = axis;
                        x = 0:.001:4*pi;
                        y = cos(x);
                        y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                        h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                        xlim([0 4*pi]);
                        title(num2str(i),'FontWeight','normal','FontSize',10);
                        if i == 1
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
                catch
                    disp('Not possible to run thetaREM Baseline vs Drug...');
                end
            end

            % DRUG 
            
            if ~isempty(ts_Drug)
                [phaseMod_Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Drug,'plotting',false,'saveMat',false);

                rippleMod = phaseMod_Drug.ripples;
                SWMod = phaseMod_Drug.SharpWave;
                thetaMod = phaseMod_Drug.theta;
                lgammaMod = phaseMod_Drug.lgamma;
                hgammaMod = phaseMod_Drug.hgamma;
                thetaRunMod = phaseMod_Drug.thetaRunMod;
                thetaREMMod = phaseMod_Drug.thetaREMMod;

                save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');

                if plt
                    % Ripples modulation
                   try
                       figure,
                       set(gcf,'Position',get(0,'ScreenSize'))
                       for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                       end 
                       saveas(gcf,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run ripples modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'SWMod');

                if plt
                    % SW Modulation
                   try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run SW modulation Baseline vs Drug');
                   end
                end

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');

                if plt
                    % Theta modulation
                   try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run theta modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');

                if plt
                    % lgamma modulation

                   try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Drug_PhaseModulation.png'])
                   catch
                       disp('Not possible to run lgamma modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');

                if plt
                    % hgamma modulation

                   try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run hgamma modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');

                if plt
                    % ThetaRun modulation
                   try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
                   end
                end

                save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');

                if plt
                    % ThetaREM modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
                    catch
                        disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
                    end
                end
            else
                rippleMod = [];
                SWMod = [];
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];
                thetaRunMod = [];
                thetaREMMod = [];
                
                save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
                
                save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
                
                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
                
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
                
                save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
                
                
                save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
                
            end
            close all;
        end
        
        % ===============================
        % 5. CELL METRICS        
        % ===============================
        
        if isempty(dir('*.cell_metrics_Baseline.cellinfo.mat')) | isempty(dir('*.cell_metrics_Drug.cellinfo.mat')) | force
            try
                if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
                    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
                    load(file.name);
                end
                    excludeManipulationIntervals = optoPulses.stimulationEpochs;
            catch
                warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
                excludeManipulationIntervals = [];
            end
            
            % Baseline
            
            cell_metrics_Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

            cell_metrics = cell_metrics_Baseline;
            save([session.general.name,'.cell_metrics_Baseline.cellinfo.mat'],'cell_metrics');
            
            % Drug
            
            if ~isempty(ts_Drug)
                cell_metrics_Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_Drug;
                save([session.general.name,'.cell_metrics_Drug.cellinfo.mat'],'cell_metrics');
            else
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_Drug.cellinfo.mat'],'cell_metrics');
            end
        end
        
        % ===============================
        % 6. ACG PEAK ( dependent upon cell_metrics);
        % ===============================
        
        minPeakTime = 15;
        
        % Baseline
        if isempty(dir('*.ACGPeak_Baseline.cellinfo.mat')) | isempty(dir('*.ACGPeak_Drug.cellinfo.mat')) | force
            
            try
                targetFile = dir('*.cell_metrics_Baseline.cellinfo.mat'); load(targetFile.name);

            catch
            end

            all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
            all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
            all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');

            pyr_color = [1 .7 .7];
            nw_color = [.7 .7 1];
            ww_color = [.7 1 1];

            spikes = loadSpikes;
            UID = spikes.UID;

            acg = cell_metrics.acg.log10;
            acg_time = cell_metrics.general.acgs.log10;
            acg_time_offset = acg_time(minPeakTime:end);
            offset = length(acg_time) - length(acg_time_offset);

            acg_smoothed = smooth(acg,10);
            acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
            acg_smoothed_norm = acg_smoothed./sum(acg_smoothed);
            acg_smoothed_offset = acg_smoothed_norm(minPeakTime:end,:);

            acgPeak = [];

            for i = 1:length(UID)
                [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
                acgPeak(i) = acg_time_offset(acgPeak_sample(i));

                [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
                acgPeak2(i) = acg_time(acgPeak_sample2(i));
            end

            acg_time_samples = acg_time;
            acg_time = log10(cell_metrics.general.acgs.log10);

            gcf = figure;
            % set(gcf,'Position',get(0,'screensize'));
            subplot(2,2,[1 2])
            hold on;
            plotFill(acg_time,acg_smoothed_norm(:,all_pyr),'Color',pyr_color);
            plotFill(acg_time,acg_smoothed_norm(:,all_nw),'Color',nw_color);
            plotFill(acg_time,acg_smoothed_norm(:,all_ww),'Color',ww_color);

            set(gca,'XTick',[(-2) (-1) 0 1])
            XTick = [-2 -1 0 1];
            XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
            set(gca,'XTickLabel',XTickLabels);
            ylabel('logACG (prob)'); xlabel('Time(s)');
            axis tight;

            subplot(2,2,3)
            hold on;
            histogram(acgPeak_sample2(all_pyr),'FaceColor',pyr_color);
            histogram(acgPeak_sample2(all_nw),'FaceColor',nw_color);
            histogram(acgPeak_sample2(all_ww),'FaceColor',ww_color);
            axis tight; ylabel('Count'); xlabel('bin number');xlim([0 60])

            subplot(2,2,4)
            hold on;
            histogram(acgPeak_sample(all_pyr)+offset,'FaceColor',pyr_color);
            histogram(acgPeak_sample(all_nw)+offset,'FaceColor',nw_color);
            histogram(acgPeak_sample(all_ww)+offset,'FaceColor',ww_color);
            axis tight; ylabel('Count'); xlabel('bin number'); xlim([0 60])

            saveas(gcf,['BaselinevsDrug\ACGPeak_Baseline.png']); 

            acgPeak = [];

            acgPeak.acg_smoothed = acg_smoothed;
            acgPeak.acg_smoothed_norm = acg_smoothed_norm;
            acgPeak.acgPeak_sample = acgPeak_sample+offset;
            acgPeak.acg_time = acg_time;
            acgPeak.acg_time_samples = acg_time_samples';

            acgPeak.acgPeak_sample2 = acgPeak_sample2;
            acgPeak.acgPeak_sample = acgPeak_sample;
            acgPeak.offset = offset;

            save([session.general.name,'.ACGPeak_Baseline.cellinfo.mat'],'acgPeak');    

            % Drug
            
            if ~isempty(ts_Drug)
                try
                    targetFile = dir('*.cell_metrics_Drug.cellinfo.mat'); load(targetFile.name);

                catch
                end

                all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
                all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
                all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');

                pyr_color = [1 .7 .7];
                nw_color = [.7 .7 1];
                ww_color = [.7 1 1];

                spikes = loadSpikes;
                UID = spikes.UID;

                acg = cell_metrics.acg.log10;
                acg_time = cell_metrics.general.acgs.log10;
                acg_time_offset = acg_time(minPeakTime:end);
                offset = length(acg_time) - length(acg_time_offset);

                acg_smoothed = smooth(acg,10);
                acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
                acg_smoothed_norm = acg_smoothed./sum(acg_smoothed);
                acg_smoothed_offset = acg_smoothed_norm(minPeakTime:end,:);

                acgPeak = [];

                for i = 1:length(UID)
                    [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
                    acgPeak(i) = acg_time_offset(acgPeak_sample(i));

                    [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
                    acgPeak2(i) = acg_time(acgPeak_sample2(i));
                end

                acg_time_samples = acg_time;
                acg_time = log10(cell_metrics.general.acgs.log10);

                gcf = figure;
                % set(gcf,'Position',get(0,'screensize'));
                subplot(2,2,[1 2])
                hold on;
                plotFill(acg_time,acg_smoothed_norm(:,all_pyr),'Color',pyr_color);
                plotFill(acg_time,acg_smoothed_norm(:,all_nw),'Color',nw_color);
                plotFill(acg_time,acg_smoothed_norm(:,all_ww),'Color',ww_color);

                set(gca,'XTick',[(-2) (-1) 0 1])
                XTick = [-2 -1 0 1];
                XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
                set(gca,'XTickLabel',XTickLabels);
                ylabel('logACG (prob)'); xlabel('Time(s)');
                axis tight;

                subplot(2,2,3)
                hold on;
                histogram(acgPeak_sample2(all_pyr),'FaceColor',pyr_color);
                histogram(acgPeak_sample2(all_nw),'FaceColor',nw_color);
                histogram(acgPeak_sample2(all_ww),'FaceColor',ww_color);
                axis tight; ylabel('Count'); xlabel('bin number');xlim([0 60])

                subplot(2,2,4)
                hold on;
                histogram(acgPeak_sample(all_pyr)+offset,'FaceColor',pyr_color);
                histogram(acgPeak_sample(all_nw)+offset,'FaceColor',nw_color);
                histogram(acgPeak_sample(all_ww)+offset,'FaceColor',ww_color);
                axis tight; ylabel('Count'); xlabel('bin number'); xlim([0 60])

                saveas(gcf,['BaselinevsDrug\ACGPeak_Drug.png']); 

                acgPeak = [];

                acgPeak.acg_smoothed = acg_smoothed;
                acgPeak.acg_smoothed_norm = acg_smoothed_norm;
                acgPeak.acgPeak_sample = acgPeak_sample+offset;
                acgPeak.acg_time = acg_time;
                acgPeak.acg_time_samples = acg_time_samples';

                acgPeak.acgPeak_sample2 = acgPeak_sample2;
                acgPeak.acgPeak_sample = acgPeak_sample;
                acgPeak.offset = offset;

                save([session.general.name,'.ACGPeak_Drug.cellinfo.mat'],'acgPeak');
            else
                
                acgPeak = [];
                
                acgPeak.acg_smoothed = [];
                acgPeak.acg_smoothed_norm = [];
                acgPeak.acgPeak_sample = [];
                acgPeak.acg_time = [];
                acgPeak.acg_time_samples = [];

                acgPeak.acgPeak_sample2 = [];
                acgPeak.acgPeak_sample = [];
                acgPeak.offset = [];
                
                save([session.general.name,'.ACGPeak_Drug.cellinfo.mat'],'acgPeak');
            end
        end
        
        
        % ===============================
        % 7. AVERAGE CCG
        % ==============================
        
        if isempty(dir('*.averageCCG_Baseline.cellinfo.mat')) | isempty(dir('*.averageCCG_Drug.cellinfo.mat')) | force
            
            winSizePlot = [-.3 .3];
            
            % Baseline
            
            averageCCG_Baseline = getAverageCCG('force',true,'includeIntervals',ts_Baseline,'savemat',false,'plotOpt',false,'saveFig',false);

            averageCCG = averageCCG_Baseline;
            save([session.general.name,'.averageCCG_Baseline.cellinfo.mat'],'averageCCG');

            t_ccg = averageCCG_Baseline.timestamps;
            allCcg = averageCCG_Baseline.allCcg;
            indCell = [1:size(allCcg,2)];

            gcf = figure,
            set(gcf,'Position',[200 -500 2500 1200]);
            for jj = 1:size(spikes.UID,2)
                % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                subplot(7,ceil(size(spikes.UID,2)/7),jj);
                cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
                imagesc(t_ccg,1:max(indCell)-1,cc)
                set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                hold on
                zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
                zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                plot(t_ccg, zmean,'k','LineWidth',2);
                xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                title(num2str(jj),'FontWeight','normal','FontSize',10);

                if jj == 1
                    ylabel('Cell');
                elseif jj == size(spikes.UID,2)
                    xlabel('Time (s)');
                else
                    set(gca,'YTick',[],'XTick',[]);
                end
            end
            saveas(gcf,['BaselineVsDrug\allCellsAverageCCG_Baseline.png']);
            
            gcf = figure;
            imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Baseline.ZmeanCCG,1)],...
                averageCCG_Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Grand CCG average','FontWeight','normal','FontSize',10);
            
            saveas(gcf,['BaselineVsDrug\grandCCGAverage_Baseline.png']);

            % Drug 
            
            if ~isempty(ts_Drug)
                averageCCG_Drug = getAverageCCG('force',true,'includeIntervals',ts_Drug,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Drug;
                save([session.general.name,'.averageCCG_Drug.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Drug.timestamps;
                allCcg = averageCCG_Drug.allCcg;
                indCell = [1:size(allCcg,2)];

                gcf = figure,
                set(gcf,'Position',[200 -500 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                    subplot(7,ceil(size(spikes.UID,2)/7),jj);
                    cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
                    imagesc(t_ccg,1:max(indCell)-1,cc)
                    set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                    hold on
                    zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
                    zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                    plot(t_ccg, zmean,'k','LineWidth',2);
                    xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('Cell');
                    elseif jj == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,['BaselineVsDrug\allCellsAverageCCG_Drug.png']);
                
                gcf = figure;
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Drug.ZmeanCCG,1)],...
                    averageCCG_Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gcf,['BaselineVsDrug\grandCCGAverage_Drug.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_Drug.cellinfo.mat'],'averageCCG');
            end
            
            
        end
        
        %% ====================================
        % AVERAGE CCG EXCLUDING RIPPLES =======
        % =====================================
        
        if isempty(dir('*averageCCGNoRipples_Baseline.cellinfo.mat')) | isempty(dir('*.averageCCG_drug.cellinfo.mat')) | force
            
            winSizePlot = [-.3 .3];
           
            % Baseline
            if ~isempty(dir('*ripples_Baseline.events.mat'))
                file = dir('*ripples_Baseline.events.mat');
                load(file.name);
            end
            ts_ripples_Baseline = [ripples.timestamps];
            
            averageCCG_Baseline = getAverageCCG('force',true,'includeIntervals',ts_Baseline,'excludeIntervals',ts_ripples_Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
            
            averageCCG = averageCCG_Baseline;
            save([session.general.name,'.averageCCGNoRipples_Baseline.cellinfo.mat'],'averageCCG');
            
            t_ccg = averageCCG_Baseline.timestamps;
            allCcg = averageCCG_Baseline.allCcg;
            indCell = [1:size(allCcg,2)];
            
            gcf = figure,
            set(gcf,'Position',[200 -500 2500 1200]);
            for jj = 1:size(spikes.UID,2)
                % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                subplot(7,ceil(size(spikes.UID,2)/7),jj);
                cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
                imagesc(t_ccg,1:max(indCell)-1,cc)
                set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                hold on
                zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
                zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                plot(t_ccg, zmean,'k','LineWidth',2);
                xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                title(num2str(jj),'FontWeight','normal','FontSize',10);

                if jj == 1
                    ylabel('Cell');
                elseif jj == size(spikes.UID,2)
                    xlabel('Time (s)');
                else
                    set(gca,'YTick',[],'XTick',[]);
                end
            end
            saveas(gcf,['BaselineVsDrug\allCellsAverageCCGNoRipples_Baseline.png']);
            
            gcf = figure;
            imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Baseline.ZmeanCCG,1)],...
                averageCCG_Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Grand CCG average','FontWeight','normal','FontSize',10);
            
            saveas(gcf,['BaselineVsDrug\grandCCGAverageNoRipples_Baseline.png']);
            
            % Drug
            
            if ~isempty(dir('*ripples_Drug.events.mat'))
                file = dir('*ripples_Drug.events.mat');
                load(file.name);
            end
            ts_ripples_Drug = [ripples.timestamps];
            
            if ~isempty(ts_Drug)
                averageCCG_Drug = getAverageCCG('force',true,'includeIntervals',ts_Drug,'excludeIntervals',ts_ripples_Drug,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Drug;
                save([session.general.name,'.averageCCGNoRipples_Drug.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Drug.timestamps;
                allCcg = averageCCG_Drug.allCcg;
                indCell = [1:size(allCcg,2)];

                gcf = figure,
                set(gcf,'Position',[200 -500 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    % fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                    subplot(7,ceil(size(spikes.UID,2)/7),jj);
                    cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
                    imagesc(t_ccg,1:max(indCell)-1,cc)
                    set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                    hold on
                    zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
                    zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                    plot(t_ccg, zmean,'k','LineWidth',2);
                    xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('Cell');
                    elseif jj == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,['BaselineVsDrug\allCellsAverageCCGNoRipples_Drug.png']);
                
                gcf = figure;
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Drug.ZmeanCCG,1)],...
                    averageCCG_Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gcf,['BaselineVsDrug\grandCCGAverage_Drug.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCGNoRipples_Drug.cellinfo.mat'],'averageCCG');
            end
            
            
        end
        
        
        
        %% ====================================
        % 8. SPEED CORR
        %% ====================================
        if isempty(dir('*.speedCorrs_Baseline.cellinfo.mat')) | isempty(dir('*.speedCorrs_Drug.cellinfo.mat')) | force
            
            % Baseline
            
            try
                speedCorr_Baseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);

                speedCorrs = speedCorr_Baseline;
                save([session.general.name,'.speedCorrs_Baseline.cellinfo.mat'],'speedCorrs');

            end
            
            % Drug
            if ~isempty(ts_Drug)
                try
                    speedCorr_Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);

                    speedCorrs = speedCorr_Drug;
                    save([session.general.name,'.speedCorrs_Drug.cellinfo.mat'],'speedCorrs');
                end
            else
                speedCorrs = [];
                    save([session.general.name,'.speedCorrs_Drug.cellinfo.mat'],'speedCorrs');
            end
        end
        
        
        %% =====================================
        % 9. SUMMARY
        %% ====================================
        
        % Baseline
        cd('BaselinevsDrug');
        if ~exist('Summary_Baseline.png','file')
            cd ..
            plotSummary_Baseline('excludePlot',{'spatialModulation'});
        else
            cd ..
        end
        % Drug
        if ~isempty(ts_Drug)
            cd('BaselinevsDrug')
            if ~exist('Summary_Drug.png','file')
                cd ..
                plotSummary_Drug('excludePlot',{'spatialModulation'});
            else 
                cd ..
            end  
        end
    end
    close all;
    clc;
end