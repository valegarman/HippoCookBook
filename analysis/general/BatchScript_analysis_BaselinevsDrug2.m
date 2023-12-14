%% BatchScript_analysis_BaselinevsDrug

clear; close all
targetProject= 'MK801Project';
cd('J:\data');
database_path = 'J:\data';
HCB_directory = what('GLUN3Project'); 

% sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
sessionsTable = readtable(['C:\Users\Jorge\Documents\GitHub\MK801Project',filesep,'indexedSessions_MK801Project_Drugs.csv']);
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
forceACGPeak = true;
forceBehavior = true;
forceACG = true;

for ii = 1:length(sessionsTable.SessionName)
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
        
        ts_PreSleep = [];
        ts_Baseline = [];
        ts_Drug = [];
        ts_MazeBaseline = [];
        ts_MazeDrug = [];
        ts_Maze1Baseline = [];
        ts_Maze1Drug = [];
        ts_Maze2Baseline = [];
        ts_Maze2Drug = [];
        ts_Maze3Baseline = [];
        ts_Maze3Drug = [];
        ts_LongSleepBaseline = [];
        ts_LongSleepDrug = [];
        ts_interMazeBaseline = [];
        ts_interMazeDrug = [];
        ts_injectionInterTrial = [];
        
        
         count = 1;
        for jj = 1:length(session.epochs)
            if strcmpi(session.epochs{jj}.behavioralParadigm, 'PreSleep')
                ts_PreSleep = [ts_PreSleep ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        for jj = 1:length(session.epochs)
            session.epochs{jj}.behavioralParadigm
                
            if strcmpi(session.epochs{jj}.behavioralParadigm, 'PreSleep') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze1Baseline')  | strcmpi(session.epochs{jj}.behavioralParadigm, 'InterMazeBaseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze2Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'Maze3Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm, 'LongSleepBaseline')
                
                ts_Baseline = [ts_Baseline ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze1Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'InterMazeDrug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze2Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'Maze3Drug') | strcmpi(session.epochs{jj}.behavioralParadigm , 'LongSleepDrug')
    
                ts_Drug = [ts_Drug; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        
        
        for jj = 1:length(session.epochs)
            if strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Baseline') 
                ts_Maze1Baseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze2Baseline')
                ts_Maze2Baseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze3Baseline')
                ts_Maze3Baseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];    
    
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Drug') 
                ts_Maze1Drug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
             elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze2Drug') 
                ts_Maze2Drug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze3Drug') 
                ts_Maze3Drug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        for jj = 1:length(session.epochs)
            if strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm,'Maze2Baseline') | strcmpi(session.epochs{jj}.behavioralParadigm,'Maze3Baseline')
                ts_MazeBaseline = [ts_MazeBaseline; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                   
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'Maze1Drug') | strcmpi(session.epochs{jj}.behavioralParadigm,'Maze2Drug') | strcmpi(session.epochs{jj}.behavioralParadigm,'Maze3Drug')
                ts_MazeDrug = [ts_MazeDrug; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        for jj = 1:length(session.epochs)
            if strcmpi(session.epochs{jj}.behavioralParadigm,'LongSleepBaseline')
                ts_LongSleepBaseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'LongSleepDrug')
                ts_LongSleepDrug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'interMazeBaseline')
                ts_interMazeBaseline = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm,'interMazeDrug')
                ts_interMazeDrug = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            end
        end
        
        for jj = 1:length(session.epochs)
            if strcmpi(session.epochs{jj}.behavioralParadigm, 'injectionInterTrial')
                ts_injectionInterTrial = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
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
%         try
%              if isempty(dir('*.OpenField_Baseline.mat')) | isempty(dir('*.OpenField_Drug.mat')) | forceBehavior
%                 performance = getSessionPerformance('includeIntervals',ts_Baseline);
% 
%                 OpenField = performance.OpenField;
%                 if isfield(performance,'YMaze')
%                     YMaze = performance.YMaze;
%                 end
%                 performance = OpenField;
%                 save([session.general.name,'.OpenField_Baseline.mat'],'performance');
% 
%                 if exist('YMaze','var')  
%                     performance = YMaze;
%                 else
%                     performance.meanSpeed = NaN;
%                     performance.distance = NaN;
%                 end
%                 save([session.general.name,'.YMaze_Baseline.mat'],'performance');
%              end
%             clear performance; clear OpenField; clear YMaze;
% 
%             % Drug
%            if isempty(dir('*.YMaze_Baseline.mat')) | isempty(dir('*.YMaze_Drug.mat')) | forceBehavior
%                performance = getSessionPerformance('includeIntervals',ts_Drug);
% 
%                 OpenField = performance.OpenField;
%                 if isfield(performance,'YMaze')
%                     YMaze = performance.YMaze;
%                 end
%                 performance = OpenField;
%                 save([session.general.name,'.OpenField_Drug.mat'],'performance');
% 
%                 if exist('YMaze','var')  
%                     performance = YMaze;
%                     
%                 else
%                     performance.meanSpeed = NaN;
%                     performance.distance = NaN;
%                 end
%                 save([session.general.name,'.YMaze_Drug.mat'],'performance');
%            end
%            
%            clear performance; clear OpenField; clear YMaze;
%         catch
%         end        
        
        %% ======================================
        % COMODULATION
        % =======================================
%         [comod] = crossFreqMod(lfp,[theta_passband(1):2:theta_passband(2)],[lgamma_passband(1):5:lgamma_passband(2)]);
%         [comod] = modIndex(lfp,[theta_passband(1):2:theta_passband(2)],[lgamma_passband(1):5:lgamma_passband(2)],1);

%         try
%             targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
%         catch
%             error('Not possible to compute ripples Baseline vs Drug.')
%         end
%         thetaChannel = thetaEpochs.channel;
%         
%         % Baseline
%         
%         lfp = getLFP(thetaChannel,'restrict',ts_Baseline);
%         [phaseamplitudemap_Baseline,ampfreqs_Baseline,phasecenters_Baseline] = phaseAmplitudeDist(lfp,theta_passband,lgamma_passband);
%         
%         % Drug
%         lfp = getLFP(thetaChannel,'restrict',ts_Drug);
%         [phaseamplitudemap_Drug,ampfreqs_Drug,phasecenters_Drug] = phaseAmplitudeDist(lfp,theta_passband,lgamma_passband);
%         
%         

        
        
        %% ======================================
        % FIRING RATE PROGRESSION
        % =======================================
        
        % GENERAL
%         try
%             targetFile = dir([session.general.name,'.cell_metrics.cellinfo.mat']); load (targetFile.name);
%         catch
%         end
%         
%         % Baseline OF
%         spikemat = bz_SpktToSpkmat(spikes,'dt',30,'units','rate','win',[ts_Maze1Baseline(1) ts_Maze1Baseline(2)]);
%         spikemat.UID = spikemat.UID';
%         spikemat.data = spikemat.data';
%         
%         figure;
%         set(gcf,'Position',[100 -100 2500 1200])
%         for jj = 1:size(spikes.UID,2)
%             fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
%             subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
%             bar(spikemat.data(jj,:));
%             
%             title(cell_metrics.putativeCellType{jj});
%             if jj == 1
%                 ylabel('FR (Hz)');
%                 set(gca,'XTick',[20 173],'XTickLabel',{'0','60'});
%             elseif jj == size(spikes.UID,2)
%                 set(gca,'XTick',[20 173],'XTickLabel',{'0','60'});
%             else
%                 set(gca,'YTick',[],'XTick',[]);
%             end
%         end
%         save(gcf,'SummaryFigures\spikeProgression.png');
%         
%         
%         wn = [ts_injectionInterTrial(1)-10*60 ts_injectionInterTrial(2)+60*60];
%         spikemat = bz_SpktToSpkmat(spikes,'dt',30,'units','rate','win',wn);
%         spikemat.UID = spikemat.UID';
%         spikemat.data = spikemat.data';
%         
%         figure;
%         set(gcf,'Position',[100 -100 2500 1200])
%         for jj = 1:size(spikes.UID,2)
%             fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
%             subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
%             bar(spikemat.data(jj,:));
%             
%             title(cell_metrics.putativeCellType{jj});
%             if jj == 1
%                 ylabel('FR (Hz)');
%                 set(gca,'XTick',[20 173],'XTickLabel',{'0','60'});
%             elseif jj == size(spikes.UID,2)
%                 set(gca,'XTick',[20 173],'XTickLabel',{'0','60'});
%             else
%                 set(gca,'YTick',[],'XTick',[]);
%             end
%         end
%         save(gcf,'SummaryFigures\spikeProgression.png');
%         
%         firingRateProgression = spikemat;
%         save([session.general.name,'.firingRateProgression.cellinfo.mat'],'firingRateProgression');
%         
%         % MAZE1BASELINE
%         
%         
%         wn = [ts_Maze1Baseline(2)-ts_Maze1Baseline(1)];
%         wn_step = wn/60/4;
%         window = [ts_Maze1Baseline(1) ts_Maze1Baseline(1)+wn_step; ts_Maze1Baseline(1)+wn_step ts_Maze1Baseline(1)+wn_step*2];
% 
%         spikemat_Baseline = bz_SpktToSpkmat(spikes,'dt',30,'units','rate','win',ts_Maze1Baseline);
%         spikemat_Baseline = bz_SpktToSpkmat(spikes,'dt',30,'units','rate','win',window);
%         spikemat.UID = spikemat.UID';
%         spikemat.data = spikemat.data';
%         figure;
%         bar(spikemat.data(1,:));
%         
%         
%         % MAZE1DRUG
        
        
        
        
        
        
        % ===============================
        % SPIKE POWER SPECTRUM
        % ===============================
        
%         params.Fs = 30000; params.fpass = [5 10]; params.tapers = [10 19]; p = 0.05; params.err = [1 p]; 
%         
%         ts = InIntervals(spikes.times{1},ts_Maze1Baseline);
%         ts = spikes.times{1}(ts);
%         
%         [S,f,R,Serr]=mtspectrumpt(ts,params);
%         figure;
%         plot(log10(S))
        
        % ===============================
        % SPIKE LFP COHERENCE
        % ==============================
        
%         try
%             targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
%         catch
%             error('Not possible to compute ripples Baseline vs Drug.')
%         end
%         thetaChannel = thetaEpochs.channel;
%         
%         params.fpass = [0 100]; params.tapers = [10 19]; params.err = [1 0.05]; params.Fs = 1250;
%         
%         E = ts_Maze1Baseline(1) + (ts_Maze1Baseline(2)-ts_Maze1Baseline(1))*rand(80,1);
%         win = [2 2];
%         lfp = getLFP(thetaChannel);
%         
%         datasp = createdatamatpt(spikes.times{1},E,win);
%         datasp1 = extractdatapt(spikes.times{1},E,1);
%         
%         datalfp = createdatamatc(double(lfp.data),E,1250,win);
%         datalfp1 = extractdatac(double(lfp.data),params.Fs,E);
%         
%         [C,phi,S12,S1,S2,f,zerosp,confC,phistd] = coherencycpt(datalfp1,datasp1,params);
%         figure
%         plot(a)
        
        % ===============================
        % RIPPLES CSD
        % ===============================
        
%         try
%             targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%         catch
%             error('Not possible to compute ripples Baseline vs Drug.')
%         end
%         try
%            targetFile = dir('*removedtsripples.mat'); load(targetFile.name); 
%         catch
%            warning('No removed intervals'); 
%            ts = [0 inf]; 
%         end
%         
%         samplingRate = 1250;
%         
%         % Baseline
%         ts_ripples_Baseline = find(InIntervals(ripples.peaks,ts_Baseline));
%         ripples_Baseline.timestamps = ripples.timestamps(ts_ripples_Baseline,:);
%         ripples_Baseline.peaks = ripples.peaks(ts_ripples_Baseline,:);
%         
%         
%         lfp = getLFP('all');
%         
%         twin = [0.1];
%         [evCsd_Baseline,lfpAvg_Baseline] = bz_eventCSD(lfp,ripples_Baseline.peaks,'channels',session.extracellular.electrodeGroups.channels{3},'twin',[twin twin],'plotLFP',true,'plotCSD',true);
%         
%         
%         % Drug
%         ts_ripples_Drug = find(InIntervals(ripples.peaks,ts_Drug));
%         ripples_Drug.timestamps = ripples.timestamps(ts_ripples_Drug,:);
%         ripples_Drug.peaks = ripples.peaks(ts_ripples_Drug,:);
%         
%         lfp = getLFP('all');
%         twin = [0.1];
%         [evCsd_Drug,lfpAvg_Drug] = bz_eventCSD(lfp,ripples_Drug.peaks,'channels',session.extracellular.electrodeGroups.channels{3},'twin',[twin twin],'plotLFP',true,'plotCSD',true);
%         
%         twin = [0.1 0.1];
%         twin = twin*samplingRate;
%         taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
%         figure;
%         subplot(1,2,1);
%         contourf(taxis,1:size(evCsd_Baseline.data,2),evCsd_Baseline.data',40,'LineColor','none');hold on;
%         colormap jet;
%         caxis([-40 40])
%         colormap jet; 
%         set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD Baseline'); 
%         plot([0 0],[1 size(evCsd_Baseline,2)],'--k');hold on;
%         
%         subplot(1,2,2)
%         contourf(taxis,1:size(evCsd_Drug.data,2),evCsd_Drug.data',40,'LineColor','none');hold on;
%         colormap jet;
%         caxis([-40 40])
%         colormap jet; 
%         set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD Drug'); 
%         plot([0 0],[1 size(evCsd_Drug,2)],'--k');hold on;
%         saveas(gcf,['SummaryFigures\CSDRipples.png']);
        
        
        
        % ================================
        % THETA CSD
        % ================================
%         try
%             targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
%         catch
%             error('Not possible to compute ripples Baseline vs Drug.')
%         end
%         
%         [a,b] = InIntervals(thetaEpochs.timestamps,ts_Maze1Baseline);
%         
%         lfp = getLFP(thetaEpochs.channel);
%         lfpT = bz_Filter(double(lfp.data),'passband',[6 12]);
%         [c,d] = findpeaks(lfpT);
%         
%         figure
%         plot(lfp.data(a))
%         hold on;
%         plot(lfpT(a))
%                
%         lfp = getLFP('all');
%         zz = InIntervals(d/1250,ts_Maze1Baseline);
%         
% %         [s,ss] = InIntervals(lfp.timestamps,ts);
% %         lfp.timestamps(s) = [];
% %         lfp.data(s,:) = [];
%         
%         [evCsd_ThetaBaseline,lfpAvg_ThetaBaseline] = bz_eventCSD(lfp,d(zz)/1250,'channels',session.extracellular.electrodeGroups.channels{3},'twin',[0.2 0.2],'plotLFP',true,'plotCSD',true);
%         
%         zz2 = InIntervals(d/1250,ts_Maze1Drug);
%         [evCsd_ThetaDrug,lfpAvg_ThetaDrug] = bz_eventCSD(lfp,d(zz2)/1250,'channels',session.extracellular.electrodeGroups.channels{3},'twin',[0.2 0.2],'plotLFP',true,'plotCSD',true);
%         
%         twin = [0.2 0.2];
%         twin = twin*samplingRate;
%         taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
%         figure;
%         subplot(1,2,1);
%         contourf(taxis,1:size(evCsd_ThetaBaseline.data,2),evCsd_ThetaBaseline.data',40,'LineColor','none');hold on;
%         colormap jet;
%         caxis([-10 10])
%         colormap jet; 
%         set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD Baseline'); 
%         plot([0 0],[1 size(evCsd_ThetaBaseline,2)],'--k');hold on;
%         
%         subplot(1,2,2)
%         contourf(taxis,1:size(evCsd_ThetaDrug.data,2),evCsd_ThetaDrug.data',40,'LineColor','none');hold on;
%         colormap jet;
%         caxis([-10 10])
%         colormap jet; 
%         set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD Drug'); 
%         plot([0 0],[1 size(evCsd_ThetaDrug,2)],'--k');hold on;
%         saveas(gcf,['SummaryFigures\CSDTheta.png']);
        
        % ===============================
        % DEEP-SUB      
        % ===============================
        
%         if isempty(dir('*.cell_metrics_Baseline.cellinfo.mat')) | isempty(dir('*.cell_metrics_Drug.cellinfo.mat')) | force
%             try
%                 if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
%                     file = dir([session.general.name,'.optogeneticPulses.events.mat']);
%                     load(file.name);
%                 end
%                     excludeManipulationIntervals = optoPulses.stimulationEpochs;
%             catch
%                 warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
%                 excludeManipulationIntervals = [];
%             end
%             
%             % DEEP-SUB
%             cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'saveMat',false);
%             deepSuperficialfromRipple = gui_DeepSuperficial(pwd);
%             save([session.general.name,'.cell_metrics.cellinfo.mat'],'cell_metrics');
%         end
        
        % ===========================
        % MONOSYNAPTIC CONNECTIONS
        % ===========================
        

%          cell_metrics_Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Baseline,'manualAdjustMonoSyn',true,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%         cell_metrics = cell_metrics_Baseline;
%         save([session.general.name,'.cell_metrics_Baseline.cellinfo.mat'],'cell_metrics');
        
    end
%     close all;
    clc;
end