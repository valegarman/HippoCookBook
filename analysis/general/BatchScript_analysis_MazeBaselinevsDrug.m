%% BatchScript_analysis_MazeBaselinevsDrug

clear; close all
targetProject= 'MK801Project';
cd('J:\data');
database_path = 'J:\data';
HCB_directory = what('GLUN3Project'); 

% sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
sessionsTable = readtable(['C:\Users\Jorge\Documents\GitHub\MK801Project',filesep,'indexedSessions_MK801Project_drugs.csv']);
forceReload = false;

win_resp = [-0.025 0.025];

ripple_passband = [120 200];
SW_passband = [2 10];
theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

notchFilter = false;

plt = true;

force = false;
forceLfp = false;
forceACGPeak = false;
forceBehavior = false;

for ii = 4:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        mkdir('BaselinevsDrug')
        close all;
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
                
        % Average CCG
        
        ts_Maze1Baseline = [];
        ts_Maze1Drug = [];
        
        for jj = 1:length(session.epochs)
            session.epochs{jj}.behavioralParadigm
            if strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Baseline')
                ts_Maze1Baseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Drug')
                ts_Maze1Drug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        if isempty(ts_Maze1Baseline)
            error('No Maze1 Baseline in this session...');
        end
        % =======================================
        % LFP
        % =======================================
        if isempty(dir('*.coherogram_Maze1Baseline.mat')) || forceLfp
            
            params.Fs = session.extracellular.srLfp; params.fpass = [1 200]; params.tapers = [3 5]; params.pad = 1;

            try
                targetFile = dir('*.thetaEpochs.states.mat'); load(targetFile.name);
            catch
                error('Not possible to load thetaEpochs...');
            end

            ts_theta = thetaEpochs.intervals;
            
            session.analysisTags.channel1 = thetaEpochs.channel;
            session.analysisTags.channel2 = thetaEpochs.channel;
            lfp1w = getLFP(session.analysisTags.channel1,'noPrompts',true);
            lfp2w = getLFP(session.analysisTags.channel2,'noPrompts',true);

            if notchFilter
                lfp1 = Notched(single(lfp1w.data),params.Fs,50);
                lfp2 = Notched(single(lfp2w.data),params.Fs,50);
            else
                lfp1 = single(lfp1w.data);
                lfp2 = single(lfp2w.data);
            end

            [coherogram,phase,S12,S1,S2,t,f] = cohgramc(lfp1,lfp2,[2 1],params);

            S12(S12==0) = NaN;
            S1(S1==0) = NaN;
            S2(S2==0) = NaN;

            S12 = log10(S12); % in Db
            S1 = log10(S1); % in Db
            S2 = log10(S2);

            S12_det = detrend(S12',2)';
            S1_det = detrend(S1',2)';
            S2_det = detrend(S2',2)';

            % MAZE1 BASELINE

            t_Maze1Baseline = t(InIntervals(t,ts_Maze1Baseline));
            coherogram_Maze1Baseline = coherogram(InIntervals(t,ts_Maze1Baseline),:);
            phase_Maze1Baseline = phase(InIntervals(t,ts_Maze1Baseline),:);
            S1_Maze1Baseline = S1_det(InIntervals(t,ts_Maze1Baseline),:);
            S2_Maze1Baseline = S2_det(InIntervals(t,ts_Maze1Baseline),:);

            t_Maze1Baseline_theta = t(InIntervals(t_Maze1Baseline,ts_theta));
            coherogram_Maze1Baseline_theta = coherogram_Maze1Baseline(InIntervals(t_Maze1Baseline,ts_theta),:);
            phase_Maze1Baseline_theta = phase_Maze1Baseline(InIntervals(t_Maze1Baseline,ts_theta),:);
            S1_Maze1Baseline_theta = S1_Maze1Baseline(InIntervals(t_Maze1Baseline,ts_theta),:);
            S2_Maze1Baseline_theta = S2_Maze1Baseline(InIntervals(t_Maze1Baseline,ts_theta),:);

            cohgram.t = t_Maze1Baseline_theta;
            cohgram.f = f;
            cohgram.coherogram = coherogram_Maze1Baseline_theta;
            cohgram.S1 = S1_Maze1Baseline_theta;
            cohgram.S2 = S2_Maze1Baseline_theta;
            cohgram.phase = phase_Maze1Baseline_theta;
            
            cohgram.NonThetaEpochs.t = t_Maze1Baseline;
            cohgram.NonThetaEpochs.f = f;
            cohgram.NonThetaEpochs.coherogram = coherogram_Maze1Baseline;
            cohgram.NonThetaEpochs.S1 = S1_Maze1Baseline;
            cohgram.NonThetaEpochs.S2 = S2_Maze1Baseline;
            cohgram.NonThetaEpochs.phase = phase_Maze1Baseline;
            
            cohgram.lfp1Channel = lfp1w.channels;
            cohgram.lfp1Region = lfp1w.region;
            cohgram.lfp2Channel = lfp2w.channels;
            cohgram.lfp2Region = lfp2w.region;

            save([session.general.name,'.coherogram_Maze1Baseline.mat'],'cohgram');

            figure('position',[200 115 1300 800])
            subplot(4,6,[1 2])
            imagesc(t_Maze1Baseline_theta,f,coherogram_Maze1Baseline_theta',[-1 1]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[7 8])
            imagesc(t_Maze1Baseline_theta,f,phase_Maze1Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[13 14])
            imagesc(t_Maze1Baseline_theta,f,S1_Maze1Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[19 20])
            imagesc(t_Maze1Baseline_theta,f,S2_Maze1Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

            subplot(4,6,[3 9 15 21])
            plotFill(f,nanmean(coherogram_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[4 10 16 22])
            plotFill(f,nanmean(phase_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[5 11 17 23])
            plotFill(f,nanmean(S1_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[6 12 18 24])
            plotFill(f,nanmean(S2_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  

            saveas(gca,['BaselineVsDrug\coherogram_Maze1Baseline.png']);    
        
        
            % MAZE1 DRUG

            t_Maze1Drug = t(InIntervals(t,ts_Maze1Drug));
            coherogram_Maze1Drug = coherogram(InIntervals(t,ts_Maze1Drug),:);
            phase_Maze1Drug = phase(InIntervals(t,ts_Maze1Drug),:);
            S1_Maze1Drug = S1_det(InIntervals(t,ts_Maze1Drug),:);
            S2_Maze1Drug = S2_det(InIntervals(t,ts_Maze1Drug),:);

            t_Maze1Drug_theta = t(InIntervals(t_Maze1Drug,ts_theta));
            coherogram_Maze1Drug_theta = coherogram_Maze1Drug(InIntervals(t_Maze1Drug,ts_theta),:);
            phase_Maze1Drug_theta = phase_Maze1Drug(InIntervals(t_Maze1Drug,ts_theta),:);
            S1_Maze1Drug_theta = S1_Maze1Drug(InIntervals(t_Maze1Drug,ts_theta),:);
            S2_Maze1Drug_theta = S2_Maze1Drug(InIntervals(t_Maze1Drug,ts_theta),:);

            cohgram = [];
            cohgram.t = t_Maze1Drug_theta;
            cohgram.f = f;
            cohgram.coherogram = coherogram_Maze1Drug_theta;
            cohgram.S1 = S1_Maze1Drug_theta;
            cohgram.S2 = S2_Maze1Drug_theta;
            cohgram.phase = phase_Maze1Drug_theta;
            
            cohgram.NonThetaEpochs.t = t_Maze1Drug;
            cohgram.NonThetaEpochs.f = f;
            cohgram.NonThetaEpochs.coherogram = coherogram_Maze1Drug;
            cohgram.NonThetaEpochs.S1 = S1_Maze1Drug;
            cohgram.NonThetaEpochs.S2 = S2_Maze1Drug;
            cohgram.NonThetaEpochs.phase = phase_Maze1Drug;
            
            cohgram.lfp1Channel = lfp1w.channels;
            cohgram.lfp1Region = lfp1w.region;
            cohgram.lfp2Channel = lfp2w.channels;
            cohgram.lfp2Region = lfp2w.region;

            save([session.general.name,'.coherogram_Maze1Drug.mat'],'cohgram');

            figure('position',[200 115 1300 800])
            subplot(4,6,[1 2])
            imagesc(t_Maze1Drug_theta,f,coherogram_Maze1Drug_theta',[-1 1]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[7 8])
            imagesc(t_Maze1Drug_theta,f,phase_Maze1Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[13 14])
            imagesc(t_Maze1Drug_theta,f,S1_Maze1Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[19 20])
            imagesc(t_Maze1Drug_theta,f,S2_Maze1Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

            subplot(4,6,[3 9 15 21])
            plotFill(f,nanmean(coherogram_Maze1Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[4 10 16 22])
            plotFill(f,nanmean(phase_Maze1Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[5 11 17 23])
            plotFill(f,nanmean(S1_Maze1Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[6 12 18 24])
            plotFill(f,nanmean(S2_Maze1Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  

            saveas(gca,['BaselineVsDrug\coherogram_Maze1Drug.png']);   
            
        end
    
    
            
        % =================================
        % 4. PHASE MODULATION
        % =================================
        
        if isempty(dir('*.theta_6-12_Maze1Baseline.PhaseLockingData.cellinfo.mat')) | isempty(dir('*.theta_6-12_Maze1Drug.PhaseLockingData.cellinfo.mat')) | force
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
            
            [phaseMod_Maze1Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze1Baseline,'plotting',false,'saveMat',false);

            rippleMod = phaseMod_Maze1Baseline.ripples;
            SWMod = phaseMod_Maze1Baseline.SharpWave;
            thetaMod = phaseMod_Maze1Baseline.theta;
            lgammaMod = phaseMod_Maze1Baseline.lgamma;
            hgammaMod = phaseMod_Maze1Baseline.hgamma;
            thetaRunMod = phaseMod_Maze1Baseline.thetaRunMod;
            thetaREMMod = phaseMod_Maze1Baseline.thetaREMMod;

            save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
            
            if plt
                
                % Ripples modulation
               try
                   figure('position',[200 115 1300 800])
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
                   saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run ripples modulation Baseline vs Drug...');
               end
            end
            
            save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
            
            if plt
                % SW Modulation
               try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run SW modulation Baseline vs Drug');
               end
            end

            save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
            
            if plt
                % Theta modulation
               try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
            if plt
                % lgamma modulation

               try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze1Baseline_PhaseModulation.png'])
               catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
            if plt
                % hgamma modulation

               try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
               end
            end

            save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
            
            % ThetaRun modulation
            if plt
               try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
               catch
                   disp('Not possible to run thetaRun modulation Baseline vs Drug...');
               end
            end
        

            save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
            
            if plt
                % ThetaREM modulation
                try
                    figure('position',[200 115 1300 800])
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
                    saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
                catch
                    disp('Not possible to run thetaREM Baseline vs Drug...');
                end
            end

            % DRUG 
            
            if ~isempty(ts_Maze1Drug)
                [phaseMod_Maze1Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze1Drug,'plotting',false,'saveMat',false);

                rippleMod = phaseMod_Maze1Drug.ripples;
                SWMod = phaseMod_Maze1Drug.SharpWave;
                thetaMod = phaseMod_Maze1Drug.theta;
                lgammaMod = phaseMod_Maze1Drug.lgamma;
                hgammaMod = phaseMod_Maze1Drug.hgamma;
                thetaRunMod = phaseMod_Maze1Drug.thetaRunMod;
                thetaREMMod = phaseMod_Maze1Drug.thetaREMMod;

                save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');

                if plt
                    % Ripples modulation
                   try
                       figure('position',[200 115 1300 800])
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
                       saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze1Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run ripples modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'SWMod');

                if plt
                    % SW Modulation
                   try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze1Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run SW modulation Baseline vs Drug');
                   end
                end

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');

                if plt
                    % Theta modulation
                   try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run theta modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');

                if plt
                    % lgamma modulation

                   try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze1Drug_PhaseModulation.png'])
                   catch
                       disp('Not possible to run lgamma modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');

                if plt
                    % hgamma modulation

                   try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze1Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run hgamma modulation Baseline vs Drug...');
                   end
                end

                save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');

                if plt
                    % ThetaRun modulation
                   try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
                   catch
                       disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
                   end
                end

                save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');

                if plt
                    % ThetaREM modulation
                    try
                        figure('position',[200 115 1300 800])
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
                        saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
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
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
                
                save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
                
                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
                
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
                
                save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
                
                
                save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
                
            end
            close all;
        end
        
        % ===============================
        % 5. CELL METRICS        
        % ===============================
        
        if isempty(dir('*.cell_metrics_Maze1Baseline.cellinfo.mat')) | isempty(dir('*.cell_metrics_Maze1Drug.cellinfo.mat')) | force
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
            
            cell_metrics_Maze1Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze1Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

            cell_metrics = cell_metrics_Maze1Baseline;
            save([session.general.name,'.cell_metrics_Maze1Baseline.cellinfo.mat'],'cell_metrics');
            
            % Drug
            
            if ~isempty(ts_Maze1Drug)
                cell_metrics_Maze1Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze1Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_Maze1Drug;
                save([session.general.name,'.cell_metrics_Maze1Drug.cellinfo.mat'],'cell_metrics');
            else
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_Maze1Drug.cellinfo.mat'],'cell_metrics');
            end
        end
        
        % ===============================
        % 6. ACG PEAK ( dependent upon cell_metrics);
        % ===============================
        
        minPeakTime = 15;
        
        % Baseline
        if isempty(dir('*.ACGPeak_Maze1Baseline.cellinfo.mat')) | isempty(dir('*.ACGPeak_Maze1Drug.cellinfo.mat')) | forceACGPeak
            
            try
                targetFile = dir('*.cell_metrics_Maze1Baseline.cellinfo.mat'); load(targetFile.name);

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
            acgPeak2 = [];
            acgPeak_sample = [];
            acgPeak_sample2 = [];

            for i = 1:length(UID)
                [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
                acgPeak(i) = acg_time_offset(acgPeak_sample(i));

                [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
                acgPeak2(i) = acg_time(acgPeak_sample2(i));
            end

            acg_time_samples = acg_time;
            acg_time = log10(cell_metrics.general.acgs.log10);

            figure('position',[200 115 1300 800])
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

            saveas(gca,['BaselinevsDrug\ACGPeak_Maze1Baseline.png']); 

            acgPeak = [];

            acgPeak.acg_smoothed = acg_smoothed;
            acgPeak.acg_smoothed_norm = acg_smoothed_norm;
            acgPeak.acgPeak_sample = acgPeak_sample+offset;
            acgPeak.acg_time = acg_time;
            acgPeak.acg_time_samples = acg_time_samples';

            acgPeak.acgPeak_sample2 = acgPeak_sample2;
            acgPeak.acgPeak_sample = acgPeak_sample;
            acgPeak.offset = offset;

            save([session.general.name,'.ACGPeak_Maze1Baseline.cellinfo.mat'],'acgPeak');    

            % Drug
            
            if ~isempty(ts_Maze1Drug)
                try
                    targetFile = dir('*.cell_metrics_Maze1Drug.cellinfo.mat'); load(targetFile.name);

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
                acgPeak2 = [];
                acgPeak_sample = [];
                acgPeak_sample2 = [];

                for i = 1:length(UID)
                    [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
                    acgPeak(i) = acg_time_offset(acgPeak_sample(i));

                    [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
                    acgPeak2(i) = acg_time(acgPeak_sample2(i));
                end

                acg_time_samples = acg_time;
                acg_time = log10(cell_metrics.general.acgs.log10);

                figure('position',[200 115 1300 800])
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

                saveas(gca,['BaselinevsDrug\ACGPeak_Maze1Drug.png']); 

                acgPeak = [];

                acgPeak.acg_smoothed = acg_smoothed;
                acgPeak.acg_smoothed_norm = acg_smoothed_norm;
                acgPeak.acgPeak_sample = acgPeak_sample+offset;
                acgPeak.acg_time = acg_time;
                acgPeak.acg_time_samples = acg_time_samples';

                acgPeak.acgPeak_sample2 = acgPeak_sample2;
                acgPeak.acgPeak_sample = acgPeak_sample;
                acgPeak.offset = offset;

                save([session.general.name,'.ACGPeak_Maze1Drug.cellinfo.mat'],'acgPeak');
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
                
                save([session.general.name,'.ACGPeak_Maze1Drug.cellinfo.mat'],'acgPeak');
            end
        end
        
        
        % ===============================
        % 7. AVERAGE CCG
        % ==============================
        
        if isempty(dir('*.averageCCG_Maze1Baseline.cellinfo.mat')) | isempty(dir('*.averageCCG_Maze1Drug.cellinfo.mat')) | force
            
            winSizePlot = [-.3 .3];
            
            % Baseline
            
            averageCCG_Maze1Baseline = getAverageCCG('force',true,'includeIntervals',ts_Maze1Baseline,'savemat',false,'plotOpt',false,'saveFig',false);

            averageCCG = averageCCG_Maze1Baseline;
            save([session.general.name,'.averageCCG_Maze1Baseline.cellinfo.mat'],'averageCCG');

            t_ccg = averageCCG_Maze1Baseline.timestamps;
            allCcg = averageCCG_Maze1Baseline.allCcg;
            indCell = [1:size(allCcg,2)];

            figure('position',[200 115 1300 800])
            for kk = 1:size(spikes.UID,2)
                % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
                subplot(7,ceil(size(spikes.UID,2)/7),kk);
                cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
                imagesc(t_ccg,1:max(indCell)-1,cc)
                set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                hold on
                zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
                zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                plot(t_ccg, zmean,'k','LineWidth',2);
                xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                title(num2str(kk),'FontWeight','normal','FontSize',10);

                if kk == 1
                    ylabel('Cell');
                elseif kk == size(spikes.UID,2)
                    xlabel('Time (s)');
                else
                    set(gca,'YTick',[],'XTick',[]);
                end
            end
            saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze1Baseline.png']);
            
            figure('position',[200 115 1300 800])
            imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze1Baseline.ZmeanCCG,1)],...
                averageCCG_Maze1Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Grand CCG average','FontWeight','normal','FontSize',10);
            
            saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze1Baseline.png']);

            % Drug 
            
            if ~isempty(ts_Maze1Drug)
                averageCCG_Maze1Drug = getAverageCCG('force',true,'includeIntervals',ts_Maze1Drug,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Maze1Drug;
                save([session.general.name,'.averageCCG_Maze1Drug.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Maze1Drug.timestamps;
                allCcg = averageCCG_Maze1Drug.allCcg;
                indCell = [1:size(allCcg,2)];

                figure('position',[200 115 1300 800])
                for kk = 1:size(spikes.UID,2)
                    % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
                    subplot(7,ceil(size(spikes.UID,2)/7),kk);
                    cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
                    imagesc(t_ccg,1:max(indCell)-1,cc)
                    set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                    hold on
                    zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
                    zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                    plot(t_ccg, zmean,'k','LineWidth',2);
                    xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                    title(num2str(kk),'FontWeight','normal','FontSize',10);

                    if kk == 1
                        ylabel('Cell');
                    elseif kk == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze1Drug.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze1Drug.ZmeanCCG,1)],...
                    averageCCG_Maze1Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze1Drug.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_Maze1Drug.cellinfo.mat'],'averageCCG');
            end
            
            
        end
        
        %% ====================================
        % AVERAGE CCG EXCLUDING RIPPLES =======
        % =====================================
        
        if isempty(dir('*averageCCGNoRipples_Maze1Baseline.cellinfo.mat')) | isempty(dir('*.averageCCG_Maze1drug.cellinfo.mat')) | force
            
            winSizePlot = [-.3 .3];
           
            % Baseline
            if ~isempty(dir('*ripples_Baseline.events.mat'))
                file = dir('*ripples_Baseline.events.mat');
                load(file.name);
            end
            ts_ripples_Maze1Baseline = [ripples.timestamps];
            
            averageCCG_Maze1Baseline = getAverageCCG('force',true,'includeIntervals',ts_Maze1Baseline,'excludeIntervals',ts_ripples_Maze1Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
            
            averageCCG = averageCCG_Maze1Baseline;
            save([session.general.name,'.averageCCGNoRipples_Maze1Baseline.cellinfo.mat'],'averageCCG');
            
            t_ccg = averageCCG_Maze1Baseline.timestamps;
            allCcg = averageCCG_Maze1Baseline.allCcg;
            indCell = [1:size(allCcg,2)];
            
            figure('position',[200 115 1300 800])
            for kk = 1:size(spikes.UID,2)
                % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
                subplot(7,ceil(size(spikes.UID,2)/7),kk);
                cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
                imagesc(t_ccg,1:max(indCell)-1,cc)
                set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                hold on
                zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
                zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                plot(t_ccg, zmean,'k','LineWidth',2);
                xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                title(num2str(kk),'FontWeight','normal','FontSize',10);

                if kk == 1
                    ylabel('Cell');
                elseif kk == size(spikes.UID,2)
                    xlabel('Time (s)');
                else
                    set(gca,'YTick',[],'XTick',[]);
                end
            end
            saveas(gca,['BaselineVsDrug\allCellsAverageCCGNoRipples_Maze1Baseline.png']);
            
            figure('position',[200 115 1300 800])
            imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze1Baseline.ZmeanCCG,1)],...
                averageCCG_Maze1Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
            title('Grand CCG average','FontWeight','normal','FontSize',10);
            
            saveas(gca,['BaselineVsDrug\grandCCGAverageNoRipples_Maze1Baseline.png']);
            
            % Drug
            
            if ~isempty(ts_Maze1Drug)
                
                if ~isempty(dir('*ripples_Drug.events.mat'))
                    file = dir('*ripples_Drug.events.mat');
                    load(file.name);
                end
                ts_ripples_Maze1Drug = [ripples.timestamps];
            
                averageCCG_Maze1Drug = getAverageCCG('force',true,'includeIntervals',ts_Maze1Drug,'excludeIntervals',ts_ripples_Maze1Drug,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Maze1Drug;
                save([session.general.name,'.averageCCGNoRipples_Maze1Drug.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Maze1Drug.timestamps;
                allCcg = averageCCG_Maze1Drug.allCcg;
                indCell = [1:size(allCcg,2)];

                figure('position',[200 115 1300 800])
                for kk = 1:size(spikes.UID,2)
                    % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
                    subplot(7,ceil(size(spikes.UID,2)/7),kk);
                    cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
                    imagesc(t_ccg,1:max(indCell)-1,cc)
                    set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
                    hold on
                    zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
                    zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
                    plot(t_ccg, zmean,'k','LineWidth',2);
                    xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
                    title(num2str(kk),'FontWeight','normal','FontSize',10);

                    if kk == 1
                        ylabel('Cell');
                    elseif kk == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gca,['BaselineVsDrug\allCellsAverageCCGNoRipples_Maze1Drug.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze1Drug.ZmeanCCG,1)],...
                    averageCCG_Maze1Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverageNoRipples_Maze1Drug.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCGNoRipples_Maze1Drug.cellinfo.mat'],'averageCCG');
            end
            
            
        end
        
        
        
        %% ====================================
        % 8. SPEED CORR
        %% ====================================
        if isempty(dir('*.speedCorrs_Maze1Baseline.cellinfo.mat')) | isempty(dir('*.speedCorrs_Maze1Drug.cellinfo.mat')) | force
            
            % Baseline
            
            try
                speedCorr_Maze1Baseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze1Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);

                speedCorrs = speedCorr_Maze1Baseline;
                save([session.general.name,'.speedCorrs_Maze1Baseline.cellinfo.mat'],'speedCorrs');

            end
            
            % Drug
            if ~isempty(ts_Maze1Drug)
                try
                    speedCorr_Maze1Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);

                    speedCorrs = speedCorr_Maze1Drug;
                    save([session.general.name,'.speedCorrs_Maze1Drug.cellinfo.mat'],'speedCorrs');
                end
            else
                speedCorrs = [];
                    save([session.general.name,'.speedCorrs_Maze1Drug.cellinfo.mat'],'speedCorrs');
            end
        end
        
        
    end
%     close all;
    clc;
end