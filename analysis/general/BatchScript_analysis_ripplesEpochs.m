%% BatchScript_analysis_ripplesPerSubsession

clear; close all
targetProject= 'SubiculumProject';
cd('D:');
database_path = 'D:\';
HCB_directory = what('SubiculumProject'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_SubiculumProject.csv']); % the variable is called allSessions
forceReload = false;

win_resp = [-0.025 0.025];

ripple_passband = [120 200];
SW_passband = [2 10];
theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
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

        for jj = 1:length(session.epochs)
            
            try
                targetFile = dir('*ripples.events.mat'); load(targetFile.name);
            catch
                error('Not possible to compute ripples Baseline vs Drug.')
            end
            
            ts = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];                
            ts_ripples = find(InIntervals(ripples.peaks,ts));
            
                
                % ================================
                % 1. Ripples 
                % ================================

                try
                    figure('position',[200 115 1300 800]);
                    plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples,'inAxis',true); 
                    mkdir('Epochs');
                    saveas(gca,['Epochs\plotRippleChannel_',session.epochs{jj}.behavioralParadigm,'.png']);
                    
                    ripples_epoch.timestamps = ripples.timestamps(ts_ripples,:);
                    
                    ripples_epoch.peaks = ripples.peaks(ts_ripples,:);
                    ripples_epoch.peakNormedPower = ripples.peakNormedPower(ts_ripples,:);
                    ripples_epoch.stdev = ripples.stdev;

                    ripples_epoch.noise.times = ripples.noise.times; 
                    ripples_epoch.noise.peaks = ripples.noise.peaks;
                    ripples_epoch.noise.peakNormedPower = ripples.noise.peakNormedPower;

                    ripples_epoch.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
                    ripples_epoch.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
                    ripples_epoch.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
                    ripples_epoch.detectorinfo.detectionparms = ripples.detectorinfo.detectionparms;
                    ripples_epoch.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
                    ripples_epoch.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;

                    ripples_epoch.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
                    ripples_epoch.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;

                    if isfield(ripples,'eventSpikingParameter')
                        ripples_epoch.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
                        ripples_epoch.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
                        ripples_epoch.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
                    end
                    
                    % Ripples stats
                    
                    ripples_epoch.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples,:);
                    ripples_epoch.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples,:);
                    ripples_epoch.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples,:);
                    ripples_epoch.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples)';
                    ripples_epoch.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples)';
                    ripples_epoch.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples,:);
                    ripples_epoch.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples,:);

                    ripples_epoch.rippleStats.data.incidence = length(ripples_epoch.peaks) / (ts(2)-ts(1));

                    ripples_epoch.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples,:),
                    ripples_epoch.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples,:);
                    ripples_epoch.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples,:);
                    ripples_epoch.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples,:);
                    ripples_epoch.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples,:);
                    ripples_epoch.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;

                    ripples_epoch.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples);
                    ripples_epoch.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
                    ripples_epoch.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
                    ripples_epoch.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
                    ripples_epoch.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples);
                    ripples_epoch.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples);
                    ripples_epoch.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples);
                    ripples_epoch.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples,:);
                    ripples_epoch.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;

                    corrBinSize = 0.01;

                    [data,t] = CCG(ripples_epoch.peaks,ones(length(ripples_epoch.peaks),1),'binSize',corrBinSize);                    
                    ripples_epoch.rippleStats.stats.acg.data = data;
                    ripples_epoch.rippleStats.stats.acg.t = t;

                    [rho,p] = corrcoef(ripples_epoch.rippleStats.data.peakAmplitude,ripples_epoch.rippleStats.data.peakFrequency);
                    ripples_epoch.rippleStats.stats.amplitudeFrequency.rho = rho;
                    ripples_epoch.rippleStats.stats.amplitudeFrequency.p = p;

                    [rho,p] = corrcoef(ripples_epoch.rippleStats.data.duration,ripples_epoch.rippleStats.data.peakFrequency);
                    ripples_epoch.rippleStats.stats.durationFrequency.rho = rho;
                    ripples_epoch.rippleStats.stats.durationFrequency.p = p;

                    [rho,p] = corrcoef(ripples_epoch.rippleStats.data.duration,ripples_epoch.rippleStats.data.peakAmplitude);
                    ripples_epoch.rippleStats.stats.durationAmplitude.rho = rho;
                    ripples_epoch.rippleStats.stats.durationAmplitude.p = p;

                    ripples = ripples_epoch;
                    save([session.general.name,'.ripples_',session.epochs{jj}.behavioralParadigm,'.events.mat'],'ripples');
                
                catch
                    warning('Not possible lo run ripplesEpochs...');
                end
                
            % ========================
            % 2. Ripples PSTH 
            % =======================
                    
            if isempty(dir(['ripples_',session.epochs{jj}.behavioralParadigm,'.events.mat']))
                try

                    psthRipples_epoch = spikesPsth([],'eventType','ripples','restrictIntervals',ts,'numRep',500,'force',true,'min_pulsesNumber',0,'saveMat',false,'savePlot',false,'rasterPlot',false,'ratePlot',false);

                    t = psthRipples_epoch.timestamps;
                    winSizePlot = [-0.5 0.5];
                    figure('position',[200 115 1300 800])
                    subplot(1,2,1)
                    imagesc([t(1) t(end)],[1 size(psthRipples_epoch.responsecurve,1)],...
                        psthRipples_epoch.responsecurveSmooth); caxis([0 10]); colormap(jet);
                    set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                    title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
                    ylabel('Cells');
                    xlabel('Time');

                    subplot(1,2,2)
                    imagesc([t(1) t(end)],[1 size(psthRipples_epoch.responsecurve,1)],...
                        psthRipples_epoch.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
                    set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
                    title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
                    ylabel('Cells');
                    xlabel('Time');
                    mkdir('Epochs');

                    saveas(gca,['Epochs\spikesPsthRate_ripples_',session.epochs{jj}.behavioralParadigm,'.png']);

                    raster = psthRipples_epoch.raster;
                    psth = rmfield(psthRipples_epoch,'raster');

                    save([session.general.name,'.ripples_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'psth','-v7.3');
                    save([session.general.name,'.ripples_',session.epochs{jj}.behavioralParadigm,'_raster.cellinfo.mat'],'raster','-v7.3');
                catch
                    warning('Not possible to run ripplesEpochs...');
                end
            end
            
            % =============================
            % 3. SPATIAL CLASSIFICATION
            % =============================
            
%             spatialModulation = computeSpatialClassification('tint',false);
%             spatialModulation_tint = computeSpatialClassification('tint',true);
            
            % ===============================
            % 4. PHASE MODULATION
            % ===============================
                        
            phaseMod_epoch = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts,'plotting',false,'saveMat',false);

            rippleMod = phaseMod_epoch.ripples;
            SWMod = phaseMod_epoch.SharpWave;
            thetaMod = phaseMod_epoch.theta;
            lgammaMod = phaseMod_epoch.lgamma;
            hgammaMod = phaseMod_epoch.hgamma;
            thetaRunMod = phaseMod_epoch.thetaRunMod;
            thetaREMMod = phaseMod_epoch.thetaREMMod;

            save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'rippleMod');

            save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'SWMod');

            save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'thetaMod');

            save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'lgammaMod');

            save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'hgammaMod');

            save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'thetaRunMod');

            save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'thetaREMMod');

            
        end
    end
    close all;
end
