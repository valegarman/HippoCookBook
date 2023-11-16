%% BatchScript_analysis_BaselinevsDrug

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

force = true;
forceLfp = true;
forceACGPeak = true;
forceBehavior = true;
forceACG = true;
forceCFC = false;

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        mkdir('BaselinevsDrug')
        close all;
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
        % CHECKING BRAIN STATES
        % The State Editor
%         TheStateEditor_temp(session.general.name);
        
        
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
        
%         for jj = 1:length(session.epochs)
%             if strcmpi(session.epochs{jj}.behavioralParadigm, 'injectionInterTrial')
%                 ts_injectionInterTrial = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
%             end
%         end
        
        ts_Baseline = [ts_Baseline(1) ts_Baseline(end)];
        if ~isempty(ts_Drug)
            ts_Drug = [ts_Drug(1) ts_Drug(end)];
        end

        
        %% ======================================
        % CROSS-FREQUENCY COUPLING
        % =======================================
%         if isempty(dir('*.CFC.mat')) || forceCFC
%             % General
%             
% %             CFC = computeCrossFrequencyCoupling();
%             CFC = computeCrossFrequencyCouplingv2();
%             
%             % Baseline
%             if ~isempty(ts_Baseline)
%                 
%                 CFC_Baseline = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Baseline,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Baseline.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Baseline.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Baseline.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Baseline.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Baseline.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Baseline_theta.png']);
%                 
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Baseline.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Baseline.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Baseline.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Baseline.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Baseline.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Baseline_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Baseline.Comodulogram;
%                 CFC.PhaseVector = CFC_Baseline.PhaseVector;
%                 CFC.AmpVector = CFC_Baseline.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Baseline.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Baseline.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Baseline.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Baseline.png']);
%                 
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Baseline.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Baseline.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Baseline.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Baseline.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Baseline.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Baseline_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Baseline_theta.png']);
%                
%             end
%             
%             % Drug
%             
%             if ~isempty(ts_Drug)
%                 
%                 CFC_Drug = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Drug,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Drug.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Drug.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Drug.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Drug.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Drug.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Drug_theta.png']);
%                 
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Drug.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Drug.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Drug.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Drug.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Drug.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Drug_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Drug.Comodulogram;
%                 CFC.PhaseVector = CFC_Drug.PhaseVector;
%                 CFC.AmpVector = CFC_Drug.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Drug.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Drug.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Drug.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Drug.png']);
%                 
%                      % theta 
%                      
%                      CFC = [];
%                     CFC.Comodulogram = CFC_Drug.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Drug.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Drug.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Drug.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Drug.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Drug_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Drug_theta.png']);
%                
%             end
%             
%             % MazeBaseline
%             
%             if ~isempty(ts_MazeBaseline)
%                 
%                 CFC_MazeBaseline = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_MazeBaseline,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_MazeBaseline.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_MazeBaseline.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_MazeBaseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_MazeBaseline.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_MazeBaseline.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_MazeBaseline.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_MazeBaseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_MazeBaseline_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_MazeBaseline.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_MazeBaseline.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_MazeBaseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_MazeBaseline.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_MazeBaseline.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_MazeBaseline.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_MazeBaseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_MazeBaseline_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_MazeBaseline.Comodulogram;
%                 CFC.PhaseVector = CFC_MazeBaseline.PhaseVector;
%                 CFC.AmpVector = CFC_MazeBaseline.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_MazeBaseline.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_MazeBaseline.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_MazeBaseline.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_MazeBaseline.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.Comodulogram = CFC_MazeBaseline.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_MazeBaseline.theta.PhaseVector;
%                     CFC.AmpVector = CFC_MazeBaseline.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_MazeBaseline.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_MazeBaseline.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_MazeBaseline_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_MazeBaseline_theta.png']);
%                
%             end
%             
%             % MazeDrug
%             
%             if ~isempty(ts_MazeDrug)
%                 
%                 CFC_MazeDrug = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_MazeDrug,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_MazeDrug.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_MazeDrug.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_MazeDrug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_MazeDrug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_MazeDrug.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_MazeDrug.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_MazeDrug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_MazeDrug_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_MazeDrug.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_MazeDrug.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_MazeDrug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_MazeDrug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_MazeDrug.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_MazeDrug.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_MazeDrug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_MazeDrug_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_MazeDrug.Comodulogram;
%                 CFC.PhaseVector = CFC_MazeDrug.PhaseVector;
%                 CFC.AmpVector = CFC_MazeDrug.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_MazeDrug.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_MazeDrug.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_MazeDrug.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_MazeDrug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_MazeDrug.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_MazeDrug.theta.PhaseVector;
%                     CFC.AmpVector = CFC_MazeDrug.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_MazeDrug.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_MazeDrug.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_MazeDrug_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_MazeDrug_theta.png']);
%                
%             end
%             
%             % Maze1Baseline
%             
%             if ~isempty(ts_Maze1Baseline)
%                 
%                 CFC_Maze1Baseline = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze1Baseline,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze1Baseline.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze1Baseline.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze1Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze1Baseline.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze1Baseline.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze1Baseline.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze1Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze1Baseline_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze1Baseline.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze1Baseline.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze1Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze1Baseline.png']);
%                 
%                     % theta
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze1Baseline.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze1Baseline.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze1Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze1Baseline_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze1Baseline.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze1Baseline.PhaseVector;
%                 CFC.AmpVector = CFC_Maze1Baseline.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze1Baseline.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze1Baseline.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze1Baseline.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze1Baseline.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze1Baseline.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze1Baseline.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze1Baseline.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze1Baseline.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze1Baseline.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze1Baseline_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze1Baseline_theta.png']);
%                
%             end
%             
%             % Maze1Drug
%             
%             if ~isempty(ts_Maze1Drug)
%                 
%                 CFC_Maze1Drug = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze1Drug,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze1Drug.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze1Drug.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze1Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze1Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze1Drug.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze1Drug.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze1Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze1Drug_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze1Drug.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze1Drug.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze1Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze1Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze1Drug.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze1Drug.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze1Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze1Drug_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze1Drug.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze1Drug.PhaseVector;
%                 CFC.AmpVector = CFC_Maze1Drug.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze1Drug.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze1Drug.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze1Drug.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze1Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze1Drug.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze1Drug.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze1Drug.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze1Drug.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze1Drug.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze1Drug_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze1Drug_theta.png']);
%                
%             end
%             
%             % Maze2Baseline
%             
%             if ~isempty(ts_Maze2Baseline)
%                 
%                 CFC_Maze2Baseline = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze2Baseline,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze2Baseline.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze2Baseline.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze2Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze2Baseline.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze2Baseline.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze2Baseline.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze2Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze2Baseline_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze2Baseline.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze2Baseline.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze2Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze2Baseline.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze2Baseline.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze2Baseline.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze2Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze2Baseline_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze2Baseline.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze2Baseline.PhaseVector;
%                 CFC.AmpVector = CFC_Maze2Baseline.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze2Baseline.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze2Baseline.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze2Baseline.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze2Baseline.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze2Baseline.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze2Baseline.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze2Baseline.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze2Baseline.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze2Baseline.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze2Baseline_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze2Baseline_theta.png']);
%                
%             end
%             
%             % Maze2Drug
%             
%             if ~isempty(ts_Maze2Drug)
%                 
%                 CFC_Maze2Drug = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze2Drug,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze2Drug.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze2Drug.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze2Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze2Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze2Drug.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze2Drug.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze2Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze2Drug_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze2Drug.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze2Drug.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze2Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze2Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze2Drug.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze2Drug.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze2Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze2Drug_theta.png']);
%                
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze2Drug.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze2Drug.PhaseVector;
%                 CFC.AmpVector = CFC_Maze2Drug.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze2Drug.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze2Drug.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze2Drug.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze2Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze2Drug.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze2Drug.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze2Drug.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze2Drug.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze2Drug.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze2Drug_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze2Drug_theta.png']);
%                
%             end
%             
%             % Maze3Baseline
%             
%             if ~isempty(ts_Maze3Baseline)
%                 
%                 CFC_Maze3Baseline = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze3Baseline,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze3Baseline.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze3Baseline.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze3Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze3Baseline.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze3Baseline.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze3Baseline.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze3Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze3Baseline_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze3Baseline.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze3Baseline.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze3Baseline.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze3Baseline.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze3Baseline.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze3Baseline.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze3Baseline_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze3Baseline_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze3Baseline.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze3Baseline.PhaseVector;
%                 CFC.AmpVector = CFC_Maze3Baseline.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze3Baseline.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze3Baseline.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze3Baseline.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze3Baseline.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze3Baseline.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze3Baseline.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze3Baseline.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze3Baseline.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze3Baseline.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze3Baseline_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze3Baseline_theta.png']);
%             end
%             
%             % MazeDrug
%             
%             if ~isempty(ts_Maze3Drug)
%                 
%                 CFC_Maze3Drug = computeCrossFrequencyCouplingv2('restrictToIntervals',ts_Maze3Drug,'saveMat',false,'saveFig',false);
%                
%                 % lGamma
%                 CFC = [];
%                 CFC.lGamma_MI = CFC_Maze3Drug.lGamma_MI;
%                 CFC.lGamma_MeanAmp = CFC_Maze3Drug.lGamma_MeanAmp;
%                 save([session.general.name,'.lgamma_CFC_Maze3Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze3Drug.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.lGamma_MI = CFC_Maze3Drug.theta.lGamma_MI;
%                     CFC.lGamma_MeanAmp = CFC_Maze3Drug.theta.lGamma_MeanAmp;
%                     save([session.general.name,'.lgamma_CFC_Maze3Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.lGamma_MeanAmp, CFC.lGamma_MeanAmp]/ sum(CFC.lGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_lGamma_Maze3Drug_theta.png']);
%                 
%                 % hGamma
%                 CFC = [];
%                 CFC.hGamma_MI = CFC_Maze3Drug.hGamma_MI;
%                 CFC.hGamma_MeanAmp = CFC_Maze3Drug.hGamma_MeanAmp;
%                 save([session.general.name,'.hgamma_CFC_Maze3Drug.mat'],'CFC');
%                
%                 figure;
%                 bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                 set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                 ylabel('lGamma Amplitude');
%                 xlabel('Theta phase (rad)');
%                 saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze3Drug.png']);
%                 
%                     % theta 
%                     
%                     CFC = [];
%                     CFC.hGamma_MI = CFC_Maze3Drug.theta.hGamma_MI;
%                     CFC.hGamma_MeanAmp = CFC_Maze3Drug.theta.hGamma_MeanAmp;
%                     save([session.general.name,'.hgamma_CFC_Maze3Drug_theta.mat'],'CFC');
% 
%                     figure;
%                     bar(10:20:720,[CFC.hGamma_MeanAmp, CFC.hGamma_MeanAmp]/ sum(CFC.hGamma_MeanAmp),'k')
%                     set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
%                     ylabel('lGamma Amplitude');
%                     xlabel('Theta phase (rad)');
%                     saveas(gcf,['BaselineVsDrug\CFC_hGamma_Maze3Drug_theta.png']);
%                 
%                 % Comodulogram
%                 
%                 CFC = [];
%                 CFC.Comodulogram = CFC_Maze3Drug.Comodulogram;
%                 CFC.PhaseVector = CFC_Maze3Drug.PhaseVector;
%                 CFC.AmpVector = CFC_Maze3Drug.AmpVector;
%                 CFC.PhaseVector_BandWidth = CFC_Maze3Drug.PhaseVector_BandWidth;
%                 CFC.AmpVector_BandWidth = CFC_Maze3Drug.AmpVector_BandWidth;
%                 
%                 save([session.general.name,'.Comodulogram_Maze3Drug.mat'],'CFC');
%                 
%                 figure,
%                 contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                 set(gca,'fontsize',14)
%                 ylabel('Amplitude Frequency (Hz)')
%                 xlabel('Phase Frequency (Hz)')
%                 colorbar
%                 title('Comodulogram whole recording')
%                 saveas(gcf,['BaselineVsDrug\Comodulogram_Maze3Drug.png']);
%                 
%                     % theta 
%                     CFC = [];
%                     CFC.Comodulogram = CFC_Maze3Drug.theta.Comodulogram;
%                     CFC.PhaseVector = CFC_Maze3Drug.theta.PhaseVector;
%                     CFC.AmpVector = CFC_Maze3Drug.theta.AmpVector;
%                     CFC.PhaseVector_BandWidth = CFC_Maze3Drug.theta.PhaseVector_BandWidth;
%                     CFC.AmpVector_BandWidth = CFC_Maze3Drug.theta.AmpVector_BandWidth;
% 
%                     save([session.general.name,'.Comodulogram_Maze3Drug_theta.mat'],'CFC');
% 
%                     figure,
%                     contourf(CFC.PhaseVector+CFC.PhaseVector_BandWidth/2,CFC.AmpVector+CFC.AmpVector_BandWidth/2,CFC.Comodulogram',30,'lines','none')
%                     set(gca,'fontsize',14)
%                     ylabel('Amplitude Frequency (Hz)')
%                     xlabel('Phase Frequency (Hz)')
%                     colorbar
%                     title('Comodulogram whole recording')
%                     saveas(gcf,['BaselineVsDrug\Comodulogram_Maze3Drug_theta.png']);
%                
%             end
%         end
        
        %% ======================================
        % FIRING RATE PROGRESSION
        % =======================================
%         wn = [ts_injectionInterTrial-10*60 ts_injectionInterTrial+60*60];
%         spikemat = bz_SpktToSpkmat(spikes,'dt',10,'units','rate','win',wn);
        
        % =======================================
        % LFP
        % =======================================
        if isempty(dir('*.coherogram_Baseline.mat')) || forceLfp
            
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
            
            S1_aux = S1;

            S12(S12==0) = NaN;
            S1(S1==0) = NaN;
            S2(S2==0) = NaN;

            S12 = log10(S12); % in Db
            S1 = log10(S1); % in Db
            S2 = log10(S2);

            S12_det = detrend(S12',2)';
            S1_det = detrend(S1',2)';
            S2_det = detrend(S2',2)';

            % Baseline 

            t_Baseline = t(InIntervals(t,ts_Baseline));
            coherogram_Baseline = coherogram(InIntervals(t,ts_Baseline),:);
            phase_Baseline = phase(InIntervals(t,ts_Baseline),:);
            S1_Baseline = S1_det(InIntervals(t,ts_Baseline),:);
            S2_Baseline = S2_det(InIntervals(t,ts_Baseline),:);
            
            S1_Baseline_aux = S1_aux(InIntervals(t,ts_Baseline),:);
            

            t_Baseline_theta = t(InIntervals(t_Baseline,ts_theta));
            coherogram_Baseline_theta = coherogram_Baseline(InIntervals(t_Baseline,ts_theta),:);
            phase_Baseline_theta = phase_Baseline(InIntervals(t_Baseline,ts_theta),:);
            S1_Baseline_theta = S1_Baseline(InIntervals(t_Baseline,ts_theta),:);
            S2_Baseline_theta = S2_Baseline(InIntervals(t_Baseline,ts_theta),:);
            
            S1_Baseline_aux_theta = S1_Baseline_aux(InIntervals(t_Baseline,ts_theta),:);
            
            

            cohgram.t = t_Baseline_theta;
            cohgram.f = f;
            cohgram.coherogram = coherogram_Baseline_theta;
            cohgram.S1 = S1_Baseline_theta;
            cohgram.S2 = S2_Baseline_theta;
            cohgram.phase = phase_Baseline_theta;
            cohgram.S1_aux = S1_Baseline_aux_theta;
            
            cohgram.NonThetaEpochs.t = t_Baseline;
            cohgram.NonThetaEpochs.f = f;
            cohgram.NonThetaEpochs.coherogram = coherogram_Baseline;
            cohgram.NonThetaEpochs.S1 = S1_Baseline;
            cohgram.NonThetaEpochs.S2 = S2_Baseline;
            cohgram.NonThetaEpochs.phase = phase_Baseline;
            cohgram.NonThetaEpochs.S1_aux = S1_Baseline_aux;
            
            cohgram.lfp1Channel = lfp1w.channels;
            cohgram.lfp1Region = lfp1w.region;
            cohgram.lfp2Channel = lfp2w.channels;
            cohgram.lfp2Region = lfp2w.region;

            save([session.general.name,'.coherogram_Baseline.mat'],'cohgram');

            figure('position',[200 115 1300 800])
            subplot(4,6,[1 2])
            imagesc(t_Baseline_theta,f,coherogram_Baseline_theta',[-1 1]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[7 8])
            imagesc(t_Baseline_theta,f,phase_Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[13 14])
            imagesc(t_Baseline_theta,f,S1_Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[19 20])
            imagesc(t_Baseline_theta,f,S2_Baseline_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

            subplot(4,6,[3 9 15 21])
            plotFill(f,nanmean(coherogram_Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[4 10 16 22])
            plotFill(f,nanmean(phase_Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[5 11 17 23])
            plotFill(f,nanmean(S1_Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Baseline [r]'); xlabel('Freq [Hz]');  
    %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[6 12 18 24])
            plotFill(f,nanmean(S2_Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Baseline [r]'); xlabel('Freq [Hz]');  

            saveas(gca,['BaselineVsDrug\coherogram_Baseline.png']);    

            % Drug

            t_Drug = t(InIntervals(t,ts_Drug));
            coherogram_Drug = coherogram(InIntervals(t,ts_Drug),:);
            phase_Drug = phase(InIntervals(t,ts_Drug),:);
            S1_Drug = S1_det(InIntervals(t,ts_Drug),:);
            S2_Drug = S2_det(InIntervals(t,ts_Drug),:);
            
            S1_Drug_aux = S1_aux(InIntervals(t,ts_Drug),:);

            t_Drug_theta = t(InIntervals(t_Drug,ts_theta));
            coherogram_Drug_theta = coherogram_Drug(InIntervals(t_Drug,ts_theta),:);
            phase_Drug_theta = phase_Drug(InIntervals(t_Drug,ts_theta),:);
            S1_Drug_theta = S1_Drug(InIntervals(t_Drug,ts_theta),:);
            S2_Drug_theta = S2_Drug(InIntervals(t_Drug,ts_theta),:);
            
            S1_Drug_theta_aux = S1_Drug_aux(InIntervals(t_Drug,ts_theta),:);

            cohgram = [];
            cohgram.t = t_Drug_theta;
            cohgram.f = f;
            cohgram.coherogram = coherogram_Drug_theta;
            cohgram.S1 = S1_Drug_theta;
            cohgram.S2 = S2_Drug_theta;
            cohgram.phase = phase_Drug_theta;
            cohgram.S1_aux = S1_Drug_theta_aux;
            
            cohgram.NonThetaEpochs.t = t_Drug;
            cohgram.NonThetaEpochs.f = f;
            cohgram.NonThetaEpochs.coherogram = coherogram_Drug;
            cohgram.NonThetaEpochs.S1 = S1_Drug;
            cohgram.NonThetaEpochs.S2 = S2_Drug;
            cohgram.NonThetaEpochs.phase = phase_Drug;
            cohgram.NoNThetaEpochs.S1_aux = S1_Drug_aux;
            
            cohgram.lfp1Channel = lfp1w.channels;
            cohgram.lfp1Region = lfp1w.region;
            cohgram.lfp2Channel = lfp2w.channels;
            cohgram.lfp2Region = lfp2w.region;

            save([session.general.name,'.coherogram_Drug.mat'],'cohgram');

            figure('position',[200 115 1300 800])
            subplot(4,6,[1 2])
            imagesc(t_Drug_theta,f,coherogram_Drug_theta',[-1 1]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[7 8])
            imagesc(t_Drug_theta,f,phase_Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[13 14])
            imagesc(t_Drug_theta,f,S1_Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[19 20])
            imagesc(t_Drug_theta,f,S2_Drug_theta',[-1.5 1.5]);
            colormap jet
            set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
            title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

            subplot(4,6,[3 9 15 21])
            plotFill(f,nanmean(coherogram_Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[4 10 16 22])
            plotFill(f,nanmean(phase_Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

            subplot(4,6,[5 11 17 23])
            plotFill(f,nanmean(S1_Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Drug [r]'); xlabel('Freq [Hz]');  
    %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

            subplot(4,6,[6 12 18 24])
            plotFill(f,nanmean(S2_Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
            ax = axis;
            fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
            fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
            fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
            ylabel('Drug [r]'); xlabel('Freq [Hz]');  

            saveas(gca,['BaselineVsDrug\coherogram_Drug.png']);    


            % MAZE1 BASELINE
            
            if ~isempty(ts_Maze1Baseline)
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
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze1Baseline.mat'],'cohgram');
            end
        
        
            % MAZE1 DRUG
            
            if ~isempty(ts_Maze1Drug)
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
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze1Drug.mat'],'cohgram');
            end
                
            % BASELINE VS DRUG
            if ~isempty(ts_Maze1Baseline)
                figure('position',[200 115 1300 800])
                subplot(1,2,1)
                plotFill(f,nanmean(S1_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]); 
                ax = axis;
                hold on;
                plotFill(f,nanmean(S1_Maze1Drug_theta),'color',[1 0 0]); 
                xlim([1 200]); ylim([-2 2]);
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1BaselineDrug [r]'); xlabel('Freq [Hz]'); 

                subplot(1,2,2)
                plotFill(f,nanmean(S1_Maze1Baseline),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]); 
                ax = axis;
                hold on;
                plotFill(f,nanmean(S1_Maze1Drug),'color',[1 0 0]); 
                xlim([1 200]); ylim([-2 2]);
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1BaselineDrug [r]'); xlabel('Freq [Hz]'); 

                saveas(gca,['BaselineVsDrug\coherogram_Maze1BaselineDrugS1.png']); 


                figure('position',[200 115 1300 800])
                subplot(1,2,1)
                plotFill(f,nanmean(S2_Maze1Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]); 
                ax = axis;
                hold on;
                plotFill(f,nanmean(S2_Maze1Drug_theta),'color',[1 0 0]); 
                xlim([1 200]); ylim([-2 2]);
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1BaselineDrug [r]'); xlabel('Freq [Hz]'); 

                subplot(1,2,2)
                plotFill(f,nanmean(S2_Maze1Baseline),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]); 
                ax = axis;
                hold on;
                plotFill(f,nanmean(S2_Maze1Drug),'color',[1 0 0]); 
                xlim([1 200]); ylim([-2 2]);
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1BaselineDrug [r]'); xlabel('Freq [Hz]');

                saveas(gca,['BaselineVsDrug\coherogram_Maze1BaselineDrugS2.png']);
            end
            
            
            % MAZE2 BASELINE
            if ~isempty(ts_Maze2Baseline)
                t_Maze2Baseline = t(InIntervals(t,ts_Maze2Baseline));
                coherogram_Maze2Baseline = coherogram(InIntervals(t,ts_Maze2Baseline),:);
                phase_Maze2Baseline = phase(InIntervals(t,ts_Maze2Baseline),:);
                S1_Maze2Baseline = S1_det(InIntervals(t,ts_Maze2Baseline),:);
                S2_Maze2Baseline = S2_det(InIntervals(t,ts_Maze2Baseline),:);

                t_Maze2Baseline_theta = t(InIntervals(t_Maze2Baseline,ts_theta));
                coherogram_Maze2Baseline_theta = coherogram_Maze2Baseline(InIntervals(t_Maze2Baseline,ts_theta),:);
                phase_Maze2Baseline_theta = phase_Maze2Baseline(InIntervals(t_Maze2Baseline,ts_theta),:);
                S1_Maze2Baseline_theta = S1_Maze2Baseline(InIntervals(t_Maze2Baseline,ts_theta),:);
                S2_Maze2Baseline_theta = S2_Maze2Baseline(InIntervals(t_Maze2Baseline,ts_theta),:);

                cohgram.t = t_Maze2Baseline_theta;
                cohgram.f = f;
                cohgram.coherogram = coherogram_Maze2Baseline_theta;
                cohgram.S1 = S1_Maze2Baseline_theta;
                cohgram.S2 = S2_Maze2Baseline_theta;
                cohgram.phase = phase_Maze2Baseline_theta;

                cohgram.NonThetaEpochs.t = t_Maze2Baseline;
                cohgram.NonThetaEpochs.f = f;
                cohgram.NonThetaEpochs.coherogram = coherogram_Maze2Baseline;
                cohgram.NonThetaEpochs.S1 = S1_Maze2Baseline;
                cohgram.NonThetaEpochs.S2 = S2_Maze2Baseline;
                cohgram.NonThetaEpochs.phase = phase_Maze2Baseline;

                cohgram.lfp1Channel = lfp1w.channels;
                cohgram.lfp1Region = lfp1w.region;
                cohgram.lfp2Channel = lfp2w.channels;
                cohgram.lfp2Region = lfp2w.region;

                save([session.general.name,'.coherogram_Maze2Baseline.mat'],'cohgram');

                figure('position',[200 115 1300 800])
                subplot(4,6,[1 2])
                imagesc(t_Maze2Baseline_theta,f,coherogram_Maze2Baseline_theta',[-1 1]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[7 8])
                imagesc(t_Maze2Baseline_theta,f,phase_Maze2Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[13 14])
                imagesc(t_Maze2Baseline_theta,f,S1_Maze2Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[19 20])
                imagesc(t_Maze2Baseline_theta,f,S2_Maze2Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

                subplot(4,6,[3 9 15 21])
                plotFill(f,nanmean(coherogram_Maze2Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[4 10 16 22])
                plotFill(f,nanmean(phase_Maze2Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[5 11 17 23])
                plotFill(f,nanmean(S1_Maze2Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[6 12 18 24])
                plotFill(f,nanmean(S2_Maze2Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  

                saveas(gca,['BaselineVsDrug\coherogram_Maze2Baseline.png']);
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze2Baseline.mat'],'cohgram');
            end
            
            
            % MAZE2 DRUG
            if ~isempty(ts_Maze2Drug)
                t_Maze2Drug = t(InIntervals(t,ts_Maze2Drug));
                coherogram_Maze2Drug = coherogram(InIntervals(t,ts_Maze2Drug),:);
                phase_Maze2Drug = phase(InIntervals(t,ts_Maze2Drug),:);
                S1_Maze2Drug = S1_det(InIntervals(t,ts_Maze2Drug),:);
                S2_Maze2Drug = S2_det(InIntervals(t,ts_Maze2Drug),:);

                t_Maze2Drug_theta = t(InIntervals(t_Maze2Drug,ts_theta));
                coherogram_Maze2Drug_theta = coherogram_Maze2Drug(InIntervals(t_Maze2Drug,ts_theta),:);
                phase_Maze2Drug_theta = phase_Maze2Drug(InIntervals(t_Maze2Drug,ts_theta),:);
                S1_Maze2Drug_theta = S1_Maze2Drug(InIntervals(t_Maze2Drug,ts_theta),:);
                S2_Maze2Drug_theta = S2_Maze2Drug(InIntervals(t_Maze2Drug,ts_theta),:);

                cohgram = [];
                cohgram.t = t_Maze2Drug_theta;
                cohgram.f = f;
                cohgram.coherogram = coherogram_Maze2Drug_theta;
                cohgram.S1 = S1_Maze2Drug_theta;
                cohgram.S2 = S2_Maze2Drug_theta;
                cohgram.phase = phase_Maze2Drug_theta;

                cohgram.NonThetaEpochs.t = t_Maze2Drug;
                cohgram.NonThetaEpochs.f = f;
                cohgram.NonThetaEpochs.coherogram = coherogram_Maze2Drug;
                cohgram.NonThetaEpochs.S1 = S1_Maze2Drug;
                cohgram.NonThetaEpochs.S2 = S2_Maze2Drug;
                cohgram.NonThetaEpochs.phase = phase_Maze2Drug;

                cohgram.lfp1Channel = lfp1w.channels;
                cohgram.lfp1Region = lfp1w.region;
                cohgram.lfp2Channel = lfp2w.channels;
                cohgram.lfp2Region = lfp2w.region;

                save([session.general.name,'.coherogram_Maze2Drug.mat'],'cohgram');

                figure('position',[200 115 1300 800])
                subplot(4,6,[1 2])
                imagesc(t_Maze2Drug_theta,f,coherogram_Maze2Drug_theta',[-1 1]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[7 8])
                imagesc(t_Maze2Drug_theta,f,phase_Maze2Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[13 14])
                imagesc(t_Maze2Drug_theta,f,S1_Maze2Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[19 20])
                imagesc(t_Maze2Drug_theta,f,S2_Maze2Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

                subplot(4,6,[3 9 15 21])
                plotFill(f,nanmean(coherogram_Maze2Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[4 10 16 22])
                plotFill(f,nanmean(phase_Maze2Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[5 11 17 23])
                plotFill(f,nanmean(S1_Maze2Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[6 12 18 24])
                plotFill(f,nanmean(S2_Maze2Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  

                saveas(gca,['BaselineVsDrug\coherogram_Maze2Drug.png']); 
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze2Drug.mat'],'cohgram');
            end
            
            % MAZE3 BASELINE
            if ~isempty(ts_Maze3Baseline)
                t_Maze3Baseline = t(InIntervals(t,ts_Maze3Baseline));
                coherogram_Maze3Baseline = coherogram(InIntervals(t,ts_Maze3Baseline),:);
                phase_Maze3Baseline = phase(InIntervals(t,ts_Maze3Baseline),:);
                S1_Maze3Baseline = S1_det(InIntervals(t,ts_Maze3Baseline),:);
                S2_Maze3Baseline = S2_det(InIntervals(t,ts_Maze3Baseline),:);

                t_Maze3Baseline_theta = t(InIntervals(t_Maze3Baseline,ts_theta));
                coherogram_Maze3Baseline_theta = coherogram_Maze3Baseline(InIntervals(t_Maze3Baseline,ts_theta),:);
                phase_Maze3Baseline_theta = phase_Maze3Baseline(InIntervals(t_Maze3Baseline,ts_theta),:);
                S1_Maze3Baseline_theta = S1_Maze3Baseline(InIntervals(t_Maze3Baseline,ts_theta),:);
                S2_Maze3Baseline_theta = S2_Maze3Baseline(InIntervals(t_Maze3Baseline,ts_theta),:);

                cohgram.t = t_Maze3Baseline_theta;
                cohgram.f = f;
                cohgram.coherogram = coherogram_Maze3Baseline_theta;
                cohgram.S1 = S1_Maze3Baseline_theta;
                cohgram.S2 = S2_Maze3Baseline_theta;
                cohgram.phase = phase_Maze3Baseline_theta;

                cohgram.NonThetaEpochs.t = t_Maze3Baseline;
                cohgram.NonThetaEpochs.f = f;
                cohgram.NonThetaEpochs.coherogram = coherogram_Maze3Baseline;
                cohgram.NonThetaEpochs.S1 = S1_Maze3Baseline;
                cohgram.NonThetaEpochs.S2 = S2_Maze3Baseline;
                cohgram.NonThetaEpochs.phase = phase_Maze3Baseline;

                cohgram.lfp1Channel = lfp1w.channels;
                cohgram.lfp1Region = lfp1w.region;
                cohgram.lfp2Channel = lfp2w.channels;
                cohgram.lfp2Region = lfp2w.region;

                save([session.general.name,'.coherogram_Maze3Baseline.mat'],'cohgram');

                figure('position',[200 115 1300 800])
                subplot(4,6,[1 2])
                imagesc(t_Maze3Baseline_theta,f,coherogram_Maze3Baseline_theta',[-1 1]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[7 8])
                imagesc(t_Maze3Baseline_theta,f,phase_Maze3Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[13 14])
                imagesc(t_Maze3Baseline_theta,f,S1_Maze3Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[19 20])
                imagesc(t_Maze3Baseline_theta,f,S2_Maze3Baseline_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

                subplot(4,6,[3 9 15 21])
                plotFill(f,nanmean(coherogram_Maze3Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[4 10 16 22])
                plotFill(f,nanmean(phase_Maze3Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[5 11 17 23])
                plotFill(f,nanmean(S1_Maze3Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  
        %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[6 12 18 24])
                plotFill(f,nanmean(S2_Maze3Baseline_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Baseline [r]'); xlabel('Freq [Hz]');  

                saveas(gca,['BaselineVsDrug\coherogram_Maze3Baseline.png']);
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze3Baseline.mat'],'cohgram');
            end
            
            
            % MAZE3 DRUG
            if ~isempty(ts_Maze3Drug)
                t_Maze3Drug = t(InIntervals(t,ts_Maze3Drug));
                coherogram_Maze3Drug = coherogram(InIntervals(t,ts_Maze3Drug),:);
                phase_Maze3Drug = phase(InIntervals(t,ts_Maze3Drug),:);
                S1_Maze3Drug = S1_det(InIntervals(t,ts_Maze3Drug),:);
                S2_Maze3Drug = S2_det(InIntervals(t,ts_Maze3Drug),:);

                t_Maze3Drug_theta = t(InIntervals(t_Maze3Drug,ts_theta));
                coherogram_Maze3Drug_theta = coherogram_Maze3Drug(InIntervals(t_Maze3Drug,ts_theta),:);
                phase_Maze3Drug_theta = phase_Maze3Drug(InIntervals(t_Maze3Drug,ts_theta),:);
                S1_Maze3Drug_theta = S1_Maze3Drug(InIntervals(t_Maze3Drug,ts_theta),:);
                S2_Maze3Drug_theta = S2_Maze3Drug(InIntervals(t_Maze3Drug,ts_theta),:);

                cohgram = [];
                cohgram.t = t_Maze3Drug_theta;
                cohgram.f = f;
                cohgram.coherogram = coherogram_Maze3Drug_theta;
                cohgram.S1 = S1_Maze3Drug_theta;
                cohgram.S2 = S2_Maze3Drug_theta;
                cohgram.phase = phase_Maze3Drug_theta;

                cohgram.NonThetaEpochs.t = t_Maze3Drug;
                cohgram.NonThetaEpochs.f = f;
                cohgram.NonThetaEpochs.coherogram = coherogram_Maze3Drug;
                cohgram.NonThetaEpochs.S1 = S1_Maze3Drug;
                cohgram.NonThetaEpochs.S2 = S2_Maze3Drug;
                cohgram.NonThetaEpochs.phase = phase_Maze3Drug;

                cohgram.lfp1Channel = lfp1w.channels;
                cohgram.lfp1Region = lfp1w.region;
                cohgram.lfp2Channel = lfp2w.channels;
                cohgram.lfp2Region = lfp2w.region;

                save([session.general.name,'.coherogram_Maze3Drug.mat'],'cohgram');

                figure('position',[200 115 1300 800])
                subplot(4,6,[1 2])
                imagesc(t_Maze3Drug_theta,f,coherogram_Maze3Drug_theta',[-1 1]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[7 8])
                imagesc(t_Maze3Drug_theta,f,phase_Maze3Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Phase Coherence Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[13 14])
                imagesc(t_Maze3Drug_theta,f,S1_Maze3Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[19 20])
                imagesc(t_Maze3Drug_theta,f,S2_Maze3Drug_theta',[-1.5 1.5]);
                colormap jet
                set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
                title(['Ch: ', num2str(session.analysisTags.channel2) , ' ' , lfp2w.region]); 

                subplot(4,6,[3 9 15 21])
                plotFill(f,nanmean(coherogram_Maze3Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[4 10 16 22])
                plotFill(f,nanmean(phase_Maze3Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-1 1]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Phase Coherence (r) Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region, ' Ch: ', num2str(session.analysisTags.channel2), ' ', lfp2w.region]); 

                subplot(4,6,[5 11 17 23])
                plotFill(f,nanmean(S1_Maze3Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  
        %         title(['Ch: ', num2str(session.analysisTags.channel1) , ' ' , lfp1w.region]); 

                subplot(4,6,[6 12 18 24])
                plotFill(f,nanmean(S2_Maze3Drug_theta),'color',[0 0 0]); xlim([1 200]); ylim([-2 2]);     
                ax = axis;
                fill([theta_passband flip(theta_passband)],[ax([3 3 4 4])],[.8 .6 .6],'EdgeColor','none','FaceAlpha',.1);
                fill([lgamma_passband flip(lgamma_passband)],[ax([3 3 4 4])],[.8 .4 .4],'EdgeColor','none','FaceAlpha',.1);
                fill([hgamma_passband flip(hgamma_passband)],[ax([3 3 4 4])],[.8 .2 .2],'EdgeColor','none','FaceAlpha',.1);
                ylabel('Maze1Drug [r]'); xlabel('Freq [Hz]');  

                saveas(gca,['BaselineVsDrug\coherogram_Maze3Drug.png']); 
            else
                cohgram = [];
                save([session.general.name,'.coherogram_Maze3Drug.mat'],'cohgram');
            end

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
        
        % =================================
        % 2. RIPPLES 
        % =================================
        
%         if isempty(dir('*.ripples_Baseline.events.mat')) | isempty(dir('*.ripples_Drug.events.mat')) | force
%             
%             % PreSleep
%             try
%                 targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%             catch
%                 error('Not possible to compute ripples Baseline vs Drug.')
%             end
%             try
%                 ts_ripples_PreSleep = find(InIntervals(ripples.peaks,ts_PreSleep));
%                 figure('position',[200 115 1300 800])
%                 plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_PreSleep,'inAxis',true);
%                 mkdir('BaselinevsDrug');
%                 saveas(gca,['BaselinevsDrug\plotRippleChannel_PreSleep.png']);
% 
% 
%                 ripples_PreSleep.timestamps = ripples.timestamps(ts_ripples_PreSleep,:);
% 
%                 ripples_PreSleep.peaks = ripples.peaks(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.peakNormedPower = ripples.peakNormedPower(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.stdev = ripples.stdev;
% 
%                 ripples_PreSleep.noise.times = ripples.noise.times; 
%                 ripples_PreSleep.noise.peaks = ripples.noise.peaks;
%                 ripples_PreSleep.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                 ripples_PreSleep.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                 ripples_PreSleep.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                 ripples_PreSleep.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                 ripples_PreSleep.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                 ripples_PreSleep.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                 ripples_PreSleep.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                 ripples_PreSleep.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                 ripples_PreSleep.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                 if isfield(ripples,'eventSpikingParameter')
%                     ripples_PreSleep.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                     ripples_PreSleep.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                     ripples_PreSleep.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                 end
% 
%                     % Ripples stats
% 
%                 ripples_PreSleep.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_PreSleep)';
%                 ripples_PreSleep.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_PreSleep)';
%                 ripples_PreSleep.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_PreSleep,:);
% 
%                 ripples_PreSleep.rippleStats.data.incidence = length(ripples_PreSleep.peaks) / (ts_PreSleep(2)-ts_PreSleep(1));
% 
%                 ripples_PreSleep.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_PreSleep,:),
%                 ripples_PreSleep.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_PreSleep);
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_PreSleep);
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_PreSleep);
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_PreSleep);
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_PreSleep,:);
%                 ripples_PreSleep.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                 corrBinSize = 0.01;
% 
%                 [data,t] = CCG(ripples_PreSleep.peaks,ones(length(ripples_PreSleep.peaks),1),'binSize',corrBinSize);                    
%                 ripples_PreSleep.rippleStats.stats.acg.data = data;
%                 ripples_PreSleep.rippleStats.stats.acg.t = t;
% 
%                 [rho,p] = corrcoef(ripples_PreSleep.rippleStats.data.peakAmplitude,ripples_PreSleep.rippleStats.data.peakFrequency);
%                 ripples_PreSleep.rippleStats.stats.amplitudeFrequency.rho = rho;
%                 ripples_PreSleep.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_PreSleep.rippleStats.data.duration,ripples_PreSleep.rippleStats.data.peakFrequency);
%                 ripples_PreSleep.rippleStats.stats.durationFrequency.rho = rho;
%                 ripples_PreSleep.rippleStats.stats.durationFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_PreSleep.rippleStats.data.duration,ripples_PreSleep.rippleStats.data.peakAmplitude);
%                 ripples_PreSleep.rippleStats.stats.durationAmplitude.rho = rho;
%                 ripples_PreSleep.rippleStats.stats.durationAmplitude.p = p;
% 
%                 ripples = ripples_PreSleep;
%                 save([session.general.name,'.ripples_PreSleep.events.mat'],'ripples');
%             catch
%                 ripples = [];
%                 save([session.general.name,'.ripples_PreSleep.events.mat'],'ripples');
%             end
%             
%             % Baseline
%             try
%                 targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%             catch
%                 error('Not possible to compute ripples Baseline vs Drug.')
%             end
%             
%             ts_ripples_Baseline = find(InIntervals(ripples.peaks,ts_Baseline));
%             figure('position',[200 115 1300 800])
%             plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_Baseline,'inAxis',true);
%             mkdir('BaselinevsDrug');
%             saveas(gca,['BaselinevsDrug\plotRippleChannel_Baseline.png']);
% 
% 
%             ripples_Baseline.timestamps = ripples.timestamps(ts_ripples_Baseline,:);
% 
%             ripples_Baseline.peaks = ripples.peaks(ts_ripples_Baseline,:);
%             ripples_Baseline.peakNormedPower = ripples.peakNormedPower(ts_ripples_Baseline,:);
%             ripples_Baseline.stdev = ripples.stdev;
% 
%             ripples_Baseline.noise.times = ripples.noise.times; 
%             ripples_Baseline.noise.peaks = ripples.noise.peaks;
%             ripples_Baseline.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%             ripples_Baseline.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%             ripples_Baseline.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%             ripples_Baseline.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%             ripples_Baseline.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%             ripples_Baseline.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%             ripples_Baseline.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%             ripples_Baseline.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%             ripples_Baseline.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
%             
%             if isfield(ripples,'eventSpikingParameter')
%                 ripples_Baseline.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                 ripples_Baseline.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                 ripples_Baseline.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%             end
% 
%                 % Ripples stats
% 
%             ripples_Baseline.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_Baseline)';
%             ripples_Baseline.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_Baseline)';
%             ripples_Baseline.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_Baseline,:);
%             
%             ripples_Baseline.rippleStats.data.incidence = length(ripples_Baseline.peaks) / (ts_Baseline(2)-ts_Baseline(1));
% 
%             ripples_Baseline.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_Baseline,:),
%             ripples_Baseline.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_Baseline);
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_Baseline);
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_Baseline);
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_Baseline);
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_Baseline,:);
%             ripples_Baseline.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%             corrBinSize = 0.01;
% 
%             [data,t] = CCG(ripples_Baseline.peaks,ones(length(ripples_Baseline.peaks),1),'binSize',corrBinSize);                    
%             ripples_Baseline.rippleStats.stats.acg.data = data;
%             ripples_Baseline.rippleStats.stats.acg.t = t;
% 
%             [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.peakAmplitude,ripples_Baseline.rippleStats.data.peakFrequency);
%             ripples_Baseline.rippleStats.stats.amplitudeFrequency.rho = rho;
%             ripples_Baseline.rippleStats.stats.amplitudeFrequency.p = p;
% 
%             [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.duration,ripples_Baseline.rippleStats.data.peakFrequency);
%             ripples_Baseline.rippleStats.stats.durationFrequency.rho = rho;
%             ripples_Baseline.rippleStats.stats.durationFrequency.p = p;
% 
%             [rho,p] = corrcoef(ripples_Baseline.rippleStats.data.duration,ripples_Baseline.rippleStats.data.peakAmplitude);
%             ripples_Baseline.rippleStats.stats.durationAmplitude.rho = rho;
%             ripples_Baseline.rippleStats.stats.durationAmplitude.p = p;
% 
%             ripples = ripples_Baseline;
%             save([session.general.name,'.ripples_Baseline.events.mat'],'ripples');
% 
% 
%             % Drug
%             if ~isempty(ts_Drug)
%                 try
%                     targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%                 catch
%                     error('Not possible to compute ripples Baseline vs Drug.')
%                 end
% 
%                 ts_ripples_Drug = find(InIntervals(ripples.peaks,ts_Drug));
% 
%                 figure('position',[200 115 1300 800])
%                 plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_Drug,'inAxis',true);
%                 mkdir('BaselinevsDrug');
%                 saveas(gca,['BaselinevsDrug\plotRippleChannel_Drug.png']);
% 
%                 ripples_Drug.timestamps = ripples.timestamps(ts_ripples_Drug,:);
% 
%                 ripples_Drug.peaks = ripples.peaks(ts_ripples_Drug,:);
%                 ripples_Drug.peakNormedPower = ripples.peakNormedPower(ts_ripples_Drug,:);
%                 ripples_Drug.stdev = ripples.stdev;
% 
%                 ripples_Drug.noise.times = ripples.noise.times; 
%                 ripples_Drug.noise.peaks = ripples.noise.peaks;
%                 ripples_Drug.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                 ripples_Drug.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                 ripples_Drug.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                 ripples_Drug.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                 ripples_Drug.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                 ripples_Drug.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                 ripples_Drug.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                 ripples_Drug.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                 ripples_Drug.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                 if isfield(ripples,'eventSpikingParameters')
%                     ripples_Drug.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                     ripples_Drug.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                     ripples_Drug.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                 end
% 
%                     % Ripples stats
% 
%                 ripples_Drug.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_Drug)';
%                 ripples_Drug.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_Drug)';
%                 ripples_Drug.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.data.incidence = length(ripples_Drug.peaks) / (ts_Drug(2)-ts_Drug(1));
% 
%                 ripples_Drug.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_Drug,:),
%                 ripples_Drug.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_Drug);
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_Drug);
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_Drug);
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_Drug);
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_Drug,:);
%                 ripples_Drug.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                 corrBinSize = 0.01;
% 
%                 [data,t] = CCG(ripples_Drug.peaks,ones(length(ripples_Drug.peaks),1),'binSize',corrBinSize);                    
%                 ripples_Drug.rippleStats.stats.acg.data = data;
%                 ripples_Drug.rippleStats.stats.acg.t = t;
% 
%                 [rho,p] = corrcoef(ripples_Drug.rippleStats.data.peakAmplitude,ripples_Drug.rippleStats.data.peakFrequency);
%                 ripples_Drug.rippleStats.stats.amplitudeFrequency.rho = rho;
%                 ripples_Drug.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_Drug.rippleStats.data.duration,ripples_Drug.rippleStats.data.peakFrequency);
%                 ripples_Drug.rippleStats.stats.durationFrequency.rho = rho;
%                 ripples_Drug.rippleStats.stats.durationFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_Drug.rippleStats.data.duration,ripples_Drug.rippleStats.data.peakAmplitude);
%                 ripples_Drug.rippleStats.stats.durationAmplitude.rho = rho;
%                 ripples_Drug.rippleStats.stats.durationAmplitude.p = p;
% 
%                 ripples = ripples_Drug;
%                 save([session.general.name,'.ripples_Drug.events.mat'],'ripples');
%             else
%                 ripples = [];
%                 save([session.general.name,'.ripples_Drug.events.mat'],'ripples');
%             end
%             
%             % LongSleepBaseline
%             if ~isempty(ts_LongSleepBaseline)
%                 try
%                     targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%                 catch
%                     error('Not possible to compute ripples Baseline vs Drug.')
%                 end
% 
%                 ts_ripples_LongSleepBaseline = find(InIntervals(ripples.peaks,ts_LongSleepBaseline));
%                 figure('position',[200 115 1300 800])
%                 plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_LongSleepBaseline,'inAxis',true);
%                 mkdir('BaselinevsDrug');
%                 saveas(gca,['BaselinevsDrug\plotRippleChannel_LongSleepBaseline.png']);
% 
% 
%                 ripples_LongSleepBaseline.timestamps = ripples.timestamps(ts_ripples_LongSleepBaseline,:);
% 
%                 ripples_LongSleepBaseline.peaks = ripples.peaks(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.peakNormedPower = ripples.peakNormedPower(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.stdev = ripples.stdev;
% 
%                 ripples_LongSleepBaseline.noise.times = ripples.noise.times; 
%                 ripples_LongSleepBaseline.noise.peaks = ripples.noise.peaks;
%                 ripples_LongSleepBaseline.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                 ripples_LongSleepBaseline.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                 ripples_LongSleepBaseline.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                 ripples_LongSleepBaseline.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                 ripples_LongSleepBaseline.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                 ripples_LongSleepBaseline.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                 ripples_LongSleepBaseline.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                 ripples_LongSleepBaseline.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                 ripples_LongSleepBaseline.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                 if isfield(ripples,'eventSpikingParameter')
%                     ripples_LongSleepBaseline.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                     ripples_LongSleepBaseline.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                     ripples_LongSleepBaseline.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                 end
% 
%                     % Ripples stats
% 
%                 ripples_LongSleepBaseline.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_LongSleepBaseline)';
%                 ripples_LongSleepBaseline.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_LongSleepBaseline)';
%                 ripples_LongSleepBaseline.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_LongSleepBaseline,:);
% 
%                 ripples_LongSleepBaseline.rippleStats.data.incidence = length(ripples_LongSleepBaseline.peaks) / (ts_LongSleepBaseline(2)-ts_LongSleepBaseline(1));
% 
%                 ripples_LongSleepBaseline.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_LongSleepBaseline,:),
%                 ripples_LongSleepBaseline.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_LongSleepBaseline);
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_LongSleepBaseline);
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_LongSleepBaseline);
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_LongSleepBaseline);
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_LongSleepBaseline,:);
%                 ripples_LongSleepBaseline.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                 corrBinSize = 0.01;
% 
%                 [data,t] = CCG(ripples_LongSleepBaseline.peaks,ones(length(ripples_LongSleepBaseline.peaks),1),'binSize',corrBinSize);                    
%                 ripples_LongSleepBaseline.rippleStats.stats.acg.data = data;
%                 ripples_LongSleepBaseline.rippleStats.stats.acg.t = t;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.peakAmplitude,ripples_LongSleepBaseline.rippleStats.data.peakFrequency);
%                 ripples_LongSleepBaseline.rippleStats.stats.amplitudeFrequency.rho = rho;
%                 ripples_LongSleepBaseline.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.duration,ripples_LongSleepBaseline.rippleStats.data.peakFrequency);
%                 ripples_LongSleepBaseline.rippleStats.stats.durationFrequency.rho = rho;
%                 ripples_LongSleepBaseline.rippleStats.stats.durationFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepBaseline.rippleStats.data.duration,ripples_LongSleepBaseline.rippleStats.data.peakAmplitude);
%                 ripples_LongSleepBaseline.rippleStats.stats.durationAmplitude.rho = rho;
%                 ripples_LongSleepBaseline.rippleStats.stats.durationAmplitude.p = p;
% 
%                 ripples = ripples_LongSleepBaseline;
%                 save([session.general.name,'.ripples_LongSleepBaseline.events.mat'],'ripples');
%             else
%                 ripples = [];
%                 save([session.general.name,'.ripples_LongSleepBaseline.events.mat'],'ripples');
%             end
%             
%             % LongSleepDrug
%             if ~isempty(ts_LongSleepDrug)
%                 try
%                     targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%                 catch
%                     error('Not possible to compute ripples Baseline vs Drug.')
%                 end
% 
%                 ts_ripples_LongSleepDrug = find(InIntervals(ripples.peaks,ts_LongSleepDrug));
%                 figure('position',[200 115 1300 800])
%                 plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_LongSleepDrug,'inAxis',true);
%                 mkdir('BaselinevsDrug');
%                 saveas(gca,['BaselinevsDrug\plotRippleChannel_LongSleepDrug.png']);
% 
% 
%                 ripples_LongSleepDrug.timestamps = ripples.timestamps(ts_ripples_LongSleepDrug,:);
% 
%                 ripples_LongSleepDrug.peaks = ripples.peaks(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.peakNormedPower = ripples.peakNormedPower(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.stdev = ripples.stdev;
% 
%                 ripples_LongSleepDrug.noise.times = ripples.noise.times; 
%                 ripples_LongSleepDrug.noise.peaks = ripples.noise.peaks;
%                 ripples_LongSleepDrug.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                 ripples_LongSleepDrug.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                 ripples_LongSleepDrug.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                 ripples_LongSleepDrug.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                 ripples_LongSleepDrug.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                 ripples_LongSleepDrug.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                 ripples_LongSleepDrug.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                 ripples_LongSleepDrug.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                 ripples_LongSleepDrug.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                 if isfield(ripples,'eventSpikingParameter')
%                     ripples_LongSleepDrug.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                     ripples_LongSleepDrug.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                     ripples_LongSleepDrug.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                 end
% 
%                     % Ripples stats
% 
%                 ripples_LongSleepDrug.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_LongSleepDrug)';
%                 ripples_LongSleepDrug.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_LongSleepDrug)';
%                 ripples_LongSleepDrug.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_LongSleepDrug,:);
% 
%                 ripples_LongSleepDrug.rippleStats.data.incidence = length(ripples_LongSleepDrug.peaks) / (ts_LongSleepDrug(2)-ts_LongSleepDrug(1));
% 
%                 ripples_LongSleepDrug.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_LongSleepDrug,:),
%                 ripples_LongSleepDrug.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_LongSleepDrug);
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_LongSleepDrug);
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_LongSleepDrug);
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_LongSleepDrug);
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_LongSleepDrug,:);
%                 ripples_LongSleepDrug.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                 corrBinSize = 0.01;
% 
%                 [data,t] = CCG(ripples_LongSleepDrug.peaks,ones(length(ripples_LongSleepDrug.peaks),1),'binSize',corrBinSize);                    
%                 ripples_LongSleepDrug.rippleStats.stats.acg.data = data;
%                 ripples_LongSleepDrug.rippleStats.stats.acg.t = t;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.peakAmplitude,ripples_LongSleepDrug.rippleStats.data.peakFrequency);
%                 ripples_LongSleepDrug.rippleStats.stats.amplitudeFrequency.rho = rho;
%                 ripples_LongSleepDrug.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.duration,ripples_LongSleepDrug.rippleStats.data.peakFrequency);
%                 ripples_LongSleepDrug.rippleStats.stats.durationFrequency.rho = rho;
%                 ripples_LongSleepDrug.rippleStats.stats.durationFrequency.p = p;
% 
%                 [rho,p] = corrcoef(ripples_LongSleepDrug.rippleStats.data.duration,ripples_LongSleepDrug.rippleStats.data.peakAmplitude);
%                 ripples_LongSleepDrug.rippleStats.stats.durationAmplitude.rho = rho;
%                 ripples_LongSleepDrug.rippleStats.stats.durationAmplitude.p = p;
% 
%                 ripples = ripples_LongSleepDrug;
%                 save([session.general.name,'.ripples_LongSleepDrug.events.mat'],'ripples');
%             else
%                 ripples = [];
%                 save([session.general.name,'.ripples_LongSleepDrug.events.mat'],'ripples');
%             end
%             
%             % interMazeBaseline
%             if ~isempty(ts_interMazeBaseline)
%                 try
%                     targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%                 catch
%                     error('Not possible to compute ripples Baseline vs Drug.')
%                 end
%                 try
% 
%                     ts_ripples_interMazeBaseline = find(InIntervals(ripples.peaks,ts_interMazeBaseline));
%                     figure('position',[200 115 1300 800])
%                     plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_interMazeBaseline,'inAxis',true);
%                     mkdir('BaselinevsDrug');
%                     saveas(gca,['BaselinevsDrug\plotRippleChannel_interMazeBaseline.png']);
% 
% 
%                     ripples_interMazeBaseline.timestamps = ripples.timestamps(ts_ripples_interMazeBaseline,:);
% 
%                     ripples_interMazeBaseline.peaks = ripples.peaks(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.peakNormedPower = ripples.peakNormedPower(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.stdev = ripples.stdev;
% 
%                     ripples_interMazeBaseline.noise.times = ripples.noise.times; 
%                     ripples_interMazeBaseline.noise.peaks = ripples.noise.peaks;
%                     ripples_interMazeBaseline.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                     ripples_interMazeBaseline.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                     ripples_interMazeBaseline.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                     ripples_interMazeBaseline.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                     ripples_interMazeBaseline.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                     ripples_interMazeBaseline.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                     ripples_interMazeBaseline.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                     ripples_interMazeBaseline.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                     ripples_interMazeBaseline.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                     if isfield(ripples,'eventSpikingParameter')
%                         ripples_interMazeBaseline.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                         ripples_interMazeBaseline.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                         ripples_interMazeBaseline.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                     end
% 
%                         % Ripples stats
% 
%                     ripples_interMazeBaseline.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_interMazeBaseline)';
%                     ripples_interMazeBaseline.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_interMazeBaseline)';
%                     ripples_interMazeBaseline.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_interMazeBaseline,:);
% 
%                     ripples_interMazeBaseline.rippleStats.data.incidence = length(ripples_interMazeBaseline.peaks) / (ts_interMazeBaseline(2)-ts_interMazeBaseline(1));
% 
%                     ripples_interMazeBaseline.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_interMazeBaseline,:),
%                     ripples_interMazeBaseline.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_interMazeBaseline);
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_interMazeBaseline);
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_interMazeBaseline);
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_interMazeBaseline);
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_interMazeBaseline,:);
%                     ripples_interMazeBaseline.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                     corrBinSize = 0.01;
% 
%                     [data,t] = CCG(ripples_interMazeBaseline.peaks,ones(length(ripples_interMazeBaseline.peaks),1),'binSize',corrBinSize);                    
%                     ripples_interMazeBaseline.rippleStats.stats.acg.data = data;
%                     ripples_interMazeBaseline.rippleStats.stats.acg.t = t;
% 
%                     [rho,p] = corrcoef(ripples_interMazeBaseline.rippleStats.data.peakAmplitude,ripples_interMazeBaseline.rippleStats.data.peakFrequency);
%                     ripples_interMazeBaseline.rippleStats.stats.amplitudeFrequency.rho = rho;
%                     ripples_interMazeBaseline.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                     [rho,p] = corrcoef(ripples_interMazeBaseline.rippleStats.data.duration,ripples_interMazeBaseline.rippleStats.data.peakFrequency);
%                     ripples_interMazeBaseline.rippleStats.stats.durationFrequency.rho = rho;
%                     ripples_interMazeBaseline.rippleStats.stats.durationFrequency.p = p;
% 
%                     [rho,p] = corrcoef(ripples_interMazeBaseline.rippleStats.data.duration,ripples_interMazeBaseline.rippleStats.data.peakAmplitude);
%                     ripples_interMazeBaseline.rippleStats.stats.durationAmplitude.rho = rho;
%                     ripples_interMazeBaseline.rippleStats.stats.durationAmplitude.p = p;
% 
%                     ripples = ripples_interMazeBaseline;
%                     save([session.general.name,'.ripples_interMazeBaseline.events.mat'],'ripples');
%                 catch
%                     ripples = [];
%                     save([session.general.name,'.ripples_interMazeBaseline.events.mat'],'ripples');
%                 end
%             else
%                 ripples = [];
%                 save([session.general.name,'.ripples_interMazeBaseline.events.mat'],'ripples');
%             end
%             
%             % interMazeDrug
%             if ~isempty(ts_interMazeDrug)
%                 try
%                     targetFile = dir('*ripples.events.mat'); load(targetFile.name);
%                 catch
%                     error('Not possible to compute ripples Baseline vs Drug.')
%                 end
%                 try
%                     ts_ripples_interMazeDrug = find(InIntervals(ripples.peaks,ts_interMazeDrug));
%                     figure('position',[200 115 1300 800])
%                     plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples,'saveFig',false,'restrictToIntervals',ts_ripples_interMazeDrug,'inAxis',true);
%                     mkdir('BaselinevsDrug');
%                     saveas(gca,['BaselinevsDrug\plotRippleChannel_interMazeDrug.png']);
% 
% 
%                     ripples_interMazeDrug.timestamps = ripples.timestamps(ts_ripples_interMazeDrug,:);
% 
%                     ripples_interMazeDrug.peaks = ripples.peaks(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.peakNormedPower = ripples.peakNormedPower(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.stdev = ripples.stdev;
% 
%                     ripples_interMazeDrug.noise.times = ripples.noise.times; 
%                     ripples_interMazeDrug.noise.peaks = ripples.noise.peaks;
%                     ripples_interMazeDrug.noise.peakNormedPower = ripples.noise.peakNormedPower;
% 
%                     ripples_interMazeDrug.detectorinfo.detectorname = ripples.detectorinfo.detectorname;
%                     ripples_interMazeDrug.detectorinfo.detectiondate = ripples.detectorinfo.detectiondate;
%                     ripples_interMazeDrug.detectorinfo.detectionintervals = ripples.detectorinfo.detectionintervals;
%                     ripples_interMazeDrug.detectorinfo.detectionsparms = ripples.detectorinfo.detectionparms;
%                     ripples_interMazeDrug.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel;
%                     ripples_interMazeDrug.detectorinfo.noisechannel = ripples.detectorinfo.noisechannel;
% 
%                     ripples_interMazeDrug.artifactsRemovalParameters.cutpoint = ripples.artifactsRemovalParameters.cutpoint;
%                     ripples_interMazeDrug.artifactsRemovalParameters.winSize = ripples.artifactsRemovalParameters.winSize;
% 
%                     if isfield(ripples,'eventSpikingParameter')
%                         ripples_interMazeDrug.eventSpikingParameters.spikingThreshold = ripples.eventSpikingParameters.spikingThreshold;
%                         ripples_interMazeDrug.eventSpikingParameters.winSize = ripples.eventSpikingParameters.winSize;
%                         ripples_interMazeDrug.eventSpikingParameters.eventSize = ripples.eventSpikingParameters.eventSize;
%                     end
% 
%                         % Ripples stats
% 
%                     ripples_interMazeDrug.rippleStats.data.peakFrequency = ripples.rippleStats.data.peakFrequency(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.data.peakAmplitude = ripples.rippleStats.data.peakAmplitude(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.data.duration = ripples.rippleStats.data.duration(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.data.spectralEntropy = ripples.rippleStats.data.spectralEntropy(:,ts_ripples_interMazeDrug)';
%                     ripples_interMazeDrug.rippleStats.data.fastRippleIndex = ripples.rippleStats.data.fastRippleIndex(:,ts_ripples_interMazeDrug)';
%                     ripples_interMazeDrug.rippleStats.data.multiTapperFreq = ripples.rippleStats.data.multiTapperFreq(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.data.interRippleFrequency = ripples.rippleStats.data.interRippleFrequency(ts_ripples_interMazeDrug,:);
% 
%                     ripples_interMazeDrug.rippleStats.data.incidence = length(ripples_interMazeDrug.peaks) / (ts_interMazeDrug(2)-ts_interMazeDrug(1));
% 
%                     ripples_interMazeDrug.rippleStats.maps.ripples_filtered = ripples.rippleStats.maps.ripples_filtered(ts_ripples_interMazeDrug,:),
%                     ripples_interMazeDrug.rippleStats.maps.ripples_raw = ripples.rippleStats.maps.ripples_raw(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.maps.frequency = ripples.rippleStats.maps.frequency(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.maps.phase = ripples.rippleStats.maps.phase(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.maps.amplitude = ripples.rippleStats.maps.amplitude(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.maps.timestamps = ripples.rippleStats.maps.timestamps;
% 
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.TimeFreq = ripples.rippleStats.maps.multitaperSpecs.TimeFreq(:,:,ts_ripples_interMazeDrug);
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.TimeFreqM = ripples.rippleStats.maps.multitaperSpecs.TimeFreqM;
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.T = ripples.rippleStats.maps.multitaperSpecs.T;
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.Freq_range = ripples.rippleStats.maps.multitaperSpecs.Freq_range;
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.Suma = ripples.rippleStats.maps.multitaperSpecs.Suma(:,ts_ripples_interMazeDrug);
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.frippindex = ripples.rippleStats.maps.multitaperSpecs.frippindex(:,ts_ripples_interMazeDrug);
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.entropyDada = ripples.rippleStats.maps.multitaperSpecs.entropyDada(:,ts_ripples_interMazeDrug);
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.freqData = ripples.rippleStats.maps.multitaperSpecs.freqData(ts_ripples_interMazeDrug,:);
%                     ripples_interMazeDrug.rippleStats.maps.multitaperSpecs.xtime = ripples.rippleStats.maps.multitaperSpecs.xtime;
% 
%                     corrBinSize = 0.01;
% 
%                     [data,t] = CCG(ripples_interMazeDrug.peaks,ones(length(ripples_interMazeDrug.peaks),1),'binSize',corrBinSize);                    
%                     ripples_interMazeDrug.rippleStats.stats.acg.data = data;
%                     ripples_interMazeDrug.rippleStats.stats.acg.t = t;
% 
%                     [rho,p] = corrcoef(ripples_interMazeDrug.rippleStats.data.peakAmplitude,ripples_interMazeDrug.rippleStats.data.peakFrequency);
%                     ripples_interMazeDrug.rippleStats.stats.amplitudeFrequency.rho = rho;
%                     ripples_interMazeDrug.rippleStats.stats.amplitudeFrequency.p = p;
% 
%                     [rho,p] = corrcoef(ripples_interMazeDrug.rippleStats.data.duration,ripples_interMazeDrug.rippleStats.data.peakFrequency);
%                     ripples_interMazeDrug.rippleStats.stats.durationFrequency.rho = rho;
%                     ripples_interMazeDrug.rippleStats.stats.durationFrequency.p = p;
% 
%                     [rho,p] = corrcoef(ripples_interMazeDrug.rippleStats.data.duration,ripples_interMazeDrug.rippleStats.data.peakAmplitude);
%                     ripples_interMazeDrug.rippleStats.stats.durationAmplitude.rho = rho;
%                     ripples_interMazeDrug.rippleStats.stats.durationAmplitude.p = p;
% 
%                     ripples = ripples_interMazeDrug;
%                     save([session.general.name,'.ripples_interMazeDrug.events.mat'],'ripples');
%                 catch
%                     ripples = [];
%                     save([session.general.name,'.ripples_interMazeDrug.events.mat'],'ripples');
%                 end
%             else
%                 ripples = [];
%                 save([session.general.name,'.ripples_LongSleepDrug.events.mat'],'ripples');
%             end
%             
%         end
               
        % =================================
        % 3. RIPPLES psth 
        % =================================
        
%         if isempty(dir('*.ripples_Baseline_psth.cellinfo.mat')) | isempty(dir('*.ripples_Drug_psth.cellinfo.mat')) | force
%             
%             ripples = rippleMasterDetector('rippleChannel',[],'SWChannel',[],'force',false,'removeOptogeneticStimulation',true);
%             
%             % PreSleep
%             if ~isempty(ts_PreSleep)
%                 psthRipples_PreSleep = spikesPsth([],'eventType','ripples','restrictIntervals',ts_PreSleep,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_PreSleep.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_PreSleep.responsecurve,1)],...
%                     psthRipples_PreSleep.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_PreSleep.responsecurve,1)],...
%                     psthRipples_PreSleep.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_PreSleep.png']);
% 
%                 raster = psthRipples_PreSleep.raster;
%                 psth = rmfield(psthRipples_PreSleep,'raster');
% 
%                 save([session.general.name,'.ripples_PreSleep_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_PreSleep_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_PreSleep_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_PreSleep_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%             % Baseline
%     
%             psthRipples_Baseline = spikesPsth([],'eventType','ripples','restrictIntervals',ts_Baseline,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%             t = psthRipples_Baseline.timestamps;
%             winSizePlot = [-0.5 0.5];
%             figure('position',[200 115 1300 800])
%             subplot(1,2,1)
%             imagesc([t(1) t(end)],[1 size(psthRipples_Baseline.responsecurve,1)],...
%                 psthRipples_Baseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
%             set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%             title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%             ylabel('Cells');
%             xlabel('Time');
% 
%             subplot(1,2,2)
%             imagesc([t(1) t(end)],[1 size(psthRipples_Baseline.responsecurve,1)],...
%                 psthRipples_Baseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%             set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%             title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%             ylabel('Cells');
%             xlabel('Time');
%             mkdir('BaselineVsDrug');
% 
%             saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_Baseline.png']);
% 
%             raster = psthRipples_Baseline.raster;
%             psth = rmfield(psthRipples_Baseline,'raster');
% 
%             save([session.general.name,'.ripples_Baseline_psth.cellinfo.mat'],'psth','-v7.3');
%             save([session.general.name,'.ripples_Baseline_raster.cellinfo.mat'],'raster','-v7.3');
% 
%             % Drug 
%             if ~isempty(ts_Drug)
%                 psthRipples_Drug = spikesPsth([],'eventType','ripples','restrictIntervals',ts_Drug,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_Drug.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_Drug.responsecurve,1)],...
%                     psthRipples_Drug.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_Drug.responsecurve,1)],...
%                     psthRipples_Drug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_Drug.png']);
% 
%                 raster = psthRipples_Drug.raster;
%                 psth = rmfield(psthRipples_Drug,'raster');
% 
%                 save([session.general.name,'.ripples_Drug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_Drug_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_Drug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_Drug_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%             
%             % LongSleepBaseline
%             if ~isempty(ts_LongSleepBaseline)
%                 
%                 psthRipples_LongSleepBaseline = spikesPsth([],'eventType','ripples','restrictIntervals',ts_LongSleepBaseline,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_LongSleepBaseline.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepBaseline.responsecurve,1)],...
%                     psthRipples_LongSleepBaseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepBaseline.responsecurve,1)],...
%                     psthRipples_LongSleepBaseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_LongSleepBaseline.png']);
% 
%                 raster = psthRipples_Baseline.raster;
%                 psth = rmfield(psthRipples_Baseline,'raster');
% 
%                 save([session.general.name,'.ripples_LongSleepBaseline_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_LongSleepBaseline_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_LongSleepBaseline_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_LongSleepBaseline_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%             % LongSleepDrug 
%             if ~isempty(ts_LongSleepDrug)
%                 
%                 psthRipples_LongSleepDrug = spikesPsth([],'eventType','ripples','restrictIntervals',ts_LongSleepDrug,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_LongSleepDrug.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepDrug.responsecurve,1)],...
%                     psthRipples_LongSleepDrug.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_LongSleepDrug.responsecurve,1)],...
%                     psthRipples_LongSleepDrug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_LongSleepDrug.png']);
% 
%                 raster = psthRipples_LongSleepDrug.raster;
%                 psth = rmfield(psthRipples_LongSleepDrug,'raster');
% 
%                 save([session.general.name,'.ripples_LongSleepDrug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_LongSleepDrug_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_LongSleepDrug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_LongSleepDrug_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%             % interMazeBaseline
%             if ~isempty(ts_interMazeBaseline)
%                 
%                 psthRipples_interMazeBaseline = spikesPsth([],'eventType','ripples','restrictIntervals',ts_interMazeBaseline,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_interMazeBaseline.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_interMazeBaseline.responsecurve,1)],...
%                     psthRipples_interMazeBaseline.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_interMazeBaseline.responsecurve,1)],...
%                     psthRipples_interMazeBaseline.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_interMazeBaseline.png']);
% 
%                 raster = psthRipples_interMazeBaseline.raster;
%                 psth = rmfield(psthRipples_interMazeBaseline,'raster');
% 
%                 save([session.general.name,'.ripples_interMazeBaseline_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_interMazeBaseline_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_interMazeBaseline_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_interMazeBaseline_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%             % interMazeDrug 
%             if ~isempty(ts_interMazeDrug)
%                 
%                 psthRipples_interMazeDrug = spikesPsth([],'eventType','ripples','restrictIntervals',ts_interMazeDrug,'numRep',500,'force',true,'minNumberOfPulses',2,'saveMat',false,'savePlot',false,'ratePlot',false);
% 
%                 t = psthRipples_LongSleepDrug.timestamps;
%                 winSizePlot = [-0.5 0.5];
%                 figure('position',[200 115 1300 800])
%                 subplot(1,2,1)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_interMazeDrug.responsecurve,1)],...
%                     psthRipples_interMazeDrug.responsecurveSmooth); caxis([0 10]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
% 
%                 subplot(1,2,2)
%                 imagesc([t(1) t(end)],[1 size(psthRipples_interMazeDrug.responsecurve,1)],...
%                     psthRipples_interMazeDrug.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
%                 ylabel('Cells');
%                 xlabel('Time');
%                 mkdir('BaselineVsDrug');
% 
%                 saveas(gca,['BaselinevsDrug\spikesPsthRate_ripples_interMazeDrug.png']);
% 
%                 raster = psthRipples_interMazeDrug.raster;
%                 psth = rmfield(psthRipples_interMazeDrug,'raster');
% 
%                 save([session.general.name,'.ripples_interMazeDrug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_interMazeDrug_raster.cellinfo.mat'],'raster','-v7.3');
%             else
%                 raster = [];
%                 psth = [];
% 
%                 save([session.general.name,'.ripples_interMazeDrug_psth.cellinfo.mat'],'psth','-v7.3');
%                 save([session.general.name,'.ripples_interMazeDrug_raster.cellinfo.mat'],'raster','-v7.3');
%             end
%             
%         end
        
        % =================================
        % 4. PHASE MODULATION
        % =================================
        
%         if isempty(dir('*.theta_6-12_Baseline.PhaseLockingData.cellinfo.mat')) | isempty(dir('*.theta_6-12_Drug.PhaseLockingData.cellinfo.mat')) | force
%             try
%                 targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
%                 rippleChannel = session.analysisTags.rippleChannel;
%                 if rippleChannel ~= ripples.detectorinfo.detectionchannel;
%                     error('session metadata and ripple channels are not the same');
%                 end
%             catch
%                 rippleChannel = ripples.detectorinfo.detectionchannel;
%             end
%             
%             try
%                 targetFile = dir('*.sharpwaves.events.mat'); load(targetFile.name);
%                 if ~isempty(targetFile)
%                     SWChannel = session.analysisTags.SWChannel;
%                     if SWChannel ~= SW.detectorinfo.detectionchannel
%                         error('session metadata and ripple channels are not the same');
%                     end
%                 else
%                     SWChannel = [];
%                 end
%             catch
%                 SWChannel = [];
%             end
%             
%             try
%                 targetFile = dir('*.thetaEpochs.states.mat'); load(targetFile.name);
%                 thetaChannel = session.analysisTags.thetaChannel;
%                 if thetaChannel ~= thetaEpochs.channel
%                     error('session metadata and theta channels are not the same');
%                 end
%             catch
%                 thetaChannel = thetaEpochs.channel;
%             end
%             
%             
%             % BASELINE
%             
%             [phaseMod_Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Baseline,'plotting',false,'saveMat',false);
% 
%             rippleMod = phaseMod_Baseline.ripples;
%             SWMod = phaseMod_Baseline.SharpWave;
%             thetaMod = phaseMod_Baseline.theta;
%             lgammaMod = phaseMod_Baseline.lgamma;
%             hgammaMod = phaseMod_Baseline.hgamma;
%             thetaRunMod = phaseMod_Baseline.thetaRunMod;
%             thetaREMMod = phaseMod_Baseline.thetaREMMod;
% 
%             save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
%             
%             if plt
%                 
%                 % Ripples modulation
%                try
%                    figure('position',[200 115 1300 800])
%                    for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                    end 
%                    saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Baseline_PhaseModulation.png']);
%                catch
%                    disp('Not possible to run ripples modulation Baseline vs Drug...');
%                end
%             end
%             
%             save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
%             
%             if plt
%                 % SW Modulation
%                try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Baseline_PhaseModulation.png']);
%                catch
%                    disp('Not possible to run SW modulation Baseline vs Drug');
%                end
%             end
% 
%             save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
%             
%             if plt
%                 % Theta modulation
%                try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
%                catch
%                    disp('Not possible to run theta modulation Baseline vs Drug...');
%                end
%             end
% 
%             save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%             if plt
%                 % lgamma modulation
% 
%                try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Baseline_PhaseModulation.png'])
%                catch
%                    disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                end
%             end
% 
%             save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%             if plt
%                 % hgamma modulation
% 
%                try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Baseline_PhaseModulation.png']);
%                catch
%                    disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                end
%             end
% 
%             save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%             
%             % ThetaRun modulation
%             if plt
%                try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
%                catch
%                    disp('Not possible to run thetaRun modulation Baseline vs Drug...');
%                end
%             end
%         
% 
%             save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                 '_Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             
%             if plt
%                 % ThetaREM modulation
%                 try
%                     figure('position',[200 115 1300 800])
%                     for i = 1:length(spikes.UID)
%                         subplot(7,ceil(size(spikes.UID,2)/7),i)
%                         area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                         hold on;
%                         ax = axis;
%                         x = 0:.001:4*pi;
%                         y = cos(x);
%                         y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                         h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                         xlim([0 4*pi]);
%                         title(num2str(i),'FontWeight','normal','FontSize',10);
%                         if i == 1
%                             ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                         elseif i == size(spikes.UID,2)
%                             set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                             xlabel('phase (rad)');
%                         else
%                             set(gca,'YTick',[],'XTick',[]);
%                         end
%                     end
%                     saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Baseline_PhaseModulation.png']);
%                 catch
%                     disp('Not possible to run thetaREM Baseline vs Drug...');
%                 end
%             end
% 
%             % DRUG 
%             
%             if ~isempty(ts_Drug)
%                 [phaseMod_Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Drug,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Drug.ripples;
%                 SWMod = phaseMod_Drug.SharpWave;
%                 thetaMod = phaseMod_Drug.theta;
%                 lgammaMod = phaseMod_Drug.lgamma;
%                 hgammaMod = phaseMod_Drug.hgamma;
%                 thetaRunMod = phaseMod_Drug.thetaRunMod;
%                 thetaREMMod = phaseMod_Drug.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
% 
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Drug_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
% 
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 if plt
%                     % ThetaRun modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Drug_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
%                     end
%                 end
%             else
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%                 
%             end
%             
%             % MAZEBASELINE
%             
%             if ~isempty(ts_MazeBaseline)
%                 
%                 [phaseMod_MazeBaseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_MazeBaseline,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_MazeBaseline.ripples;
%                 SWMod = phaseMod_MazeBaseline.SharpWave;
%                 thetaMod = phaseMod_MazeBaseline.theta;
%                 lgammaMod = phaseMod_MazeBaseline.lgamma;
%                 hgammaMod = phaseMod_MazeBaseline.hgamma;
%                 thetaRunMod = phaseMod_MazeBaseline.thetaRunMod;
%                 thetaREMMod = phaseMod_MazeBaseline.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
% 
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_MazeBaseline_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 % ThetaRun modulation
%                 if plt
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug...');
%                    end
%                 end
% 
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeBaseline_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM Baseline vs Drug...');
%                     end
%                 end
%             else
%                 
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeBaseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             % MAZEDRUG 
%             
%             if ~isempty(ts_MazeDrug)
%                 [phaseMod_MazeDrug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_MazeDrug,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_MazeDrug.ripples;
%                 SWMod = phaseMod_MazeDrug.SharpWave;
%                 thetaMod = phaseMod_MazeDrug.theta;
%                 lgammaMod = phaseMod_MazeDrug.lgamma;
%                 hgammaMod = phaseMod_MazeDrug.hgamma;
%                 thetaRunMod = phaseMod_MazeDrug.thetaRunMod;
%                 thetaREMMod = phaseMod_MazeDrug.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
% 
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_MazeDrug_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
% 
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 if plt
%                     % ThetaRun modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_MazeDrug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_MazeDrug_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
%                     end
%                 end
%             else
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             % MAZE1BASELINE
%             
%             if ~isempty(ts_Maze1Baseline)
%                 
%                 [phaseMod_Maze1Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze1Baseline,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze1Baseline.ripples;
%                 SWMod = phaseMod_Maze1Baseline.SharpWave;
%                 thetaMod = phaseMod_Maze1Baseline.theta;
%                 lgammaMod = phaseMod_Maze1Baseline.lgamma;
%                 hgammaMod = phaseMod_Maze1Baseline.hgamma;
%                 thetaRunMod = phaseMod_Maze1Baseline.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze1Baseline.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
% 
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze1Baseline_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 % ThetaRun modulation
%                 if plt
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug...');
%                    end
%                 end
% 
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Baseline_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM Baseline vs Drug...');
%                     end
%                 end
%             else
%                 
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             % MAZE1DRUG 
%             
%             if ~isempty(ts_Maze1Drug)
%                 [phaseMod_Maze1Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze1Drug,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze1Drug.ripples;
%                 SWMod = phaseMod_Maze1Drug.SharpWave;
%                 thetaMod = phaseMod_Maze1Drug.theta;
%                 lgammaMod = phaseMod_Maze1Drug.lgamma;
%                 hgammaMod = phaseMod_Maze1Drug.hgamma;
%                 thetaRunMod = phaseMod_Maze1Drug.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze1Drug.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
% 
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze1Drug_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
% 
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 if plt
%                     % ThetaRun modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze1Drug_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
%                     end
%                 end
%             else
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze1Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             % MAZE2BASELINE
%             
%             if ~isempty(ts_Maze2Baseline)
%                 
%                 [phaseMod_Maze2Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze2Baseline,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze2Baseline.ripples;
%                 SWMod = phaseMod_Maze2Baseline.SharpWave;
%                 thetaMod = phaseMod_Maze2Baseline.theta;
%                 lgammaMod = phaseMod_Maze2Baseline.lgamma;
%                 hgammaMod = phaseMod_Maze2Baseline.hgamma;
%                 thetaRunMod = phaseMod_Maze2Baseline.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze2Baseline.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
% 
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze2Baseline_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 % ThetaRun modulation
%                 if plt
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug...');
%                    end
%                 end
% 
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Baseline_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM Baseline vs Drug...');
%                     end
%                 end
%             else
%                 
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             % MAZE2DRUG 
%             
%             if ~isempty(ts_Maze2Drug)
%                 [phaseMod_Maze2Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze2Drug,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze2Drug.ripples;
%                 SWMod = phaseMod_Maze2Drug.SharpWave;
%                 thetaMod = phaseMod_Maze2Drug.theta;
%                 lgammaMod = phaseMod_Maze2Drug.lgamma;
%                 hgammaMod = phaseMod_Maze2Drug.hgamma;
%                 thetaRunMod = phaseMod_Maze2Drug.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze2Drug.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
% 
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze2Drug_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
% 
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 if plt
%                     % ThetaRun modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze2Drug_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
%                     end
%                 end
%             else
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze2Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             % MAZE3BASELINE
%             
%             if ~isempty(ts_Maze3Baseline)
%                 
%                 [phaseMod_Maze3Baseline] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze3Baseline,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze3Baseline.ripples;
%                 SWMod = phaseMod_Maze3Baseline.SharpWave;
%                 thetaMod = phaseMod_Maze3Baseline.theta;
%                 lgammaMod = phaseMod_Maze3Baseline.lgamma;
%                 hgammaMod = phaseMod_Maze3Baseline.hgamma;
%                 thetaRunMod = phaseMod_Maze3Baseline.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze3Baseline.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
% 
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze3Baseline_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 % ThetaRun modulation
%                 if plt
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug...');
%                    end
%                 end
% 
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Baseline_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM Baseline vs Drug...');
%                     end
%                 end
%             else
%                 
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Baseline.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             % MAZE3DRUG 
%             
%             if ~isempty(ts_Maze3Drug)
%                 [phaseMod_Maze3Drug] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Maze3Drug,'plotting',false,'saveMat',false);
% 
%                 rippleMod = phaseMod_Maze3Drug.ripples;
%                 SWMod = phaseMod_Maze3Drug.SharpWave;
%                 thetaMod = phaseMod_Maze3Drug.theta;
%                 lgammaMod = phaseMod_Maze3Drug.lgamma;
%                 hgammaMod = phaseMod_Maze3Drug.hgamma;
%                 thetaRunMod = phaseMod_Maze3Drug.thetaRunMod;
%                 thetaREMMod = phaseMod_Maze3Drug.thetaREMMod;
% 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
% 
%                 if plt
%                     % Ripples modulation
%                    try
%                        figure('position',[200 115 1300 800])
%                        for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                        end 
%                        saveas(gca,['BaselinevsDrug\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run ripples modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
% 
%                 if plt
%                     % SW Modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run SW modulation Baseline vs Drug');
%                    end
%                 end
% 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
% 
%                 if plt
%                     % Theta modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run theta modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
% 
%                 if plt
%                     % lgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Maze3Drug_PhaseModulation.png'])
%                    catch
%                        disp('Not possible to run lgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
% 
%                 if plt
%                     % hgamma modulation
% 
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run hgamma modulation Baseline vs Drug...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
% 
%                 if plt
%                     % ThetaRun modulation
%                    try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                    catch
%                        disp('Not possible to run thetaRun modulation Baseline vs Drug ...');
%                    end
%                 end
% 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
% 
%                 if plt
%                     % ThetaREM modulation
%                     try
%                         figure('position',[200 115 1300 800])
%                         for i = 1:length(spikes.UID)
%                             subplot(7,ceil(size(spikes.UID,2)/7),i)
%                             area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
%                             hold on;
%                             ax = axis;
%                             x = 0:.001:4*pi;
%                             y = cos(x);
%                             y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
%                             h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
%                             xlim([0 4*pi]);
%                             title(num2str(i),'FontWeight','normal','FontSize',10);
%                             if i == 1
%                                 ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
%                             elseif i == size(spikes.UID,2)
%                                 set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
%                                 xlabel('phase (rad)');
%                             else
%                                 set(gca,'YTick',[],'XTick',[]);
%                             end
%                         end
%                         saveas(gca,['BaselinevsDrug\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Maze3Drug_PhaseModulation.png']);
%                     catch
%                         disp('Not possible to run thetaREM modulation Baseline vs Drug ...');
%                     end
%                 end
%             else
%                 rippleMod = [];
%                 SWMod = [];
%                 thetaMod = [];
%                 lgammaMod = [];
%                 hgammaMod = [];
%                 thetaRunMod = [];
%                 thetaREMMod = [];
%                 
%                 save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'rippleMod');
%                 
%                 save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'SWMod');
%                 
%                 save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaMod');
%                 
%                 save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'lgammaMod');
%                 
%                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'hgammaMod');
%                 
%                 save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
%                 
%                 save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
%                     '_Maze3Drug.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
%             end
%             
%             
%             close all;
%         end
        
        % ===============================
        % 5. CELL METRICS        
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
%             % Baseline
%             
%             cell_metrics_Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%             cell_metrics = cell_metrics_Baseline;
%             save([session.general.name,'.cell_metrics_Baseline.cellinfo.mat'],'cell_metrics');
%             
%             % Drug
%             if ~isempty(ts_Drug)
%                 cell_metrics_Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Drug;
%                 save([session.general.name,'.cell_metrics_Drug.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Drug.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % MazeBaseline
%             if ~isempty(ts_MazeBaseline)
%                 cell_metrics_MazeBaseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_MazeBaseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_MazeBaseline;
%                 save([session.general.name,'.cell_metrics_MazeBaseline.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_MazeBaseline.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % MazeDrug
%             if ~isempty(ts_MazeDrug)
%                 cell_metrics_MazeDrug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_MazeDrug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_MazeDrug;
%                 save([session.general.name,'.cell_metrics_MazeDrug.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_MazeDrug.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % Maze1Baseline
%             if ~isempty(ts_Maze1Baseline)
%                 cell_metrics_Maze1Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze1Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze1Baseline;
%                 save([session.general.name,'.cell_metrics_Maze1Baseline.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze1Baseline.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % Maze1Drug
%             if ~isempty(ts_Maze1Drug)
%                 cell_metrics_Maze1Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze1Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze1Drug;
%                 save([session.general.name,'.cell_metrics_Maze1Drug.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze1Drug.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % Maze2Baseline
%             if ~isempty(ts_Maze2Baseline)
%                 cell_metrics_Maze2Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze2Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze2Baseline;
%                 save([session.general.name,'.cell_metrics_Maze2Baseline.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze2Baseline.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % Maze2Drug
%             if ~isempty(ts_Maze2Drug)
%                 cell_metrics_Maze2Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze2Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze2Drug;
%                 save([session.general.name,'.cell_metrics_Maze2Drug.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze2Drug.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % Maze3Baseline
%             if ~isempty(ts_Maze3Baseline)
%                 cell_metrics_Maze3Baseline = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze3Baseline,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze3Baseline;
%                 save([session.general.name,'.cell_metrics_Maze3Baseline.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze3Baseline.cellinfo.mat'],'cell_metrics');
%             end
%             
%             % MazeDrug
%             if ~isempty(ts_Maze3Drug)
%                 cell_metrics_Maze3Drug = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Maze3Drug,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%                 cell_metrics = cell_metrics_Maze3Drug;
%                 save([session.general.name,'.cell_metrics_Maze3Drug.cellinfo.mat'],'cell_metrics');
%             else
%                 cell_metrics = [];
%                 save([session.general.name,'.cell_metrics_Maze3Drug.cellinfo.mat'],'cell_metrics');
%             end
%             
%         end
%         
%         % ===============================
%         % 6. ACG PEAK ( dependent upon cell_metrics);
%         % ===============================
%         
%         minPeakTime = 15;
%         
%         % Baseline
%         if isempty(dir('*.ACGPeak_Baseline.cellinfo.mat')) | isempty(dir('*.ACGPeak_Drug.cellinfo.mat')) | forceACGPeak
%             
%             try
%                 targetFile = dir('*.cell_metrics_Baseline.cellinfo.mat'); load(targetFile.name);
% 
%             catch
%             end
% 
%             all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
%             all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
%             all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');
% 
%             pyr_color = [1 .7 .7];
%             nw_color = [.7 .7 1];
%             ww_color = [.7 1 1];
% 
%             spikes = loadSpikes;
%             UID = spikes.UID;
% 
%             acg = cell_metrics.acg.log10;
%             acg_time = cell_metrics.general.acgs.log10;
%             acg_time_offset = acg_time(minPeakTime:end);
%             offset = length(acg_time) - length(acg_time_offset);
% 
%             acg_smoothed = smooth(acg,10);
%             acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
%             acg_smoothed_norm = acg_smoothed./sum(acg_smoothed);
%             acg_smoothed_offset = acg_smoothed_norm(minPeakTime:end,:);
% 
%             acgPeak = [];
%             acgPeak2 = [];
%             acgPeak_sample = [];
%             acgPeak_sample2 = [];
% 
%             for i = 1:length(UID)
%                 [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
%                 acgPeak(i) = acg_time_offset(acgPeak_sample(i));
% 
%                 [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
%                 acgPeak2(i) = acg_time(acgPeak_sample2(i));
%             end
% 
%             acg_time_samples = acg_time;
%             acg_time = log10(cell_metrics.general.acgs.log10);
% 
%             figure('position',[200 115 1300 800])
%             % set(gcf,'Position',get(0,'screensize'));
%             subplot(2,2,[1 2])
%             hold on;
%             plotFill(acg_time,acg_smoothed_norm(:,all_pyr),'Color',pyr_color);
%             plotFill(acg_time,acg_smoothed_norm(:,all_nw),'Color',nw_color);
%             plotFill(acg_time,acg_smoothed_norm(:,all_ww),'Color',ww_color);
% 
%             set(gca,'XTick',[(-2) (-1) 0 1])
%             XTick = [-2 -1 0 1];
%             XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
%             set(gca,'XTickLabel',XTickLabels);
%             ylabel('logACG (prob)'); xlabel('Time(s)');
%             axis tight;
% 
%             subplot(2,2,3)
%             hold on;
%             histogram(acgPeak_sample2(all_pyr),'FaceColor',pyr_color);
%             histogram(acgPeak_sample2(all_nw),'FaceColor',nw_color);
%             histogram(acgPeak_sample2(all_ww),'FaceColor',ww_color);
%             axis tight; ylabel('Count'); xlabel('bin number');xlim([0 60])
% 
%             subplot(2,2,4)
%             hold on;
%             histogram(acgPeak_sample(all_pyr)+offset,'FaceColor',pyr_color);
%             histogram(acgPeak_sample(all_nw)+offset,'FaceColor',nw_color);
%             histogram(acgPeak_sample(all_ww)+offset,'FaceColor',ww_color);
%             axis tight; ylabel('Count'); xlabel('bin number'); xlim([0 60])
% 
%             saveas(gca,['BaselinevsDrug\ACGPeak_Baseline.png']); 
% 
%             acgPeak = [];
% 
%             acgPeak.acg_smoothed = acg_smoothed;
%             acgPeak.acg_smoothed_norm = acg_smoothed_norm;
%             acgPeak.acgPeak_sample = acgPeak_sample+offset;
%             acgPeak.acg_time = acg_time;
%             acgPeak.acg_time_samples = acg_time_samples';
% 
%             acgPeak.acgPeak_sample2 = acgPeak_sample2;
%             acgPeak.acgPeak_sample = acgPeak_sample;
%             acgPeak.offset = offset;
% 
%             save([session.general.name,'.ACGPeak_Baseline.cellinfo.mat'],'acgPeak');    
% 
%             % Drug
%             
%             if ~isempty(ts_Drug)
%                 try
%                     targetFile = dir('*.cell_metrics_Drug.cellinfo.mat'); load(targetFile.name);
% 
%                 catch
%                 end
% 
%                 all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
%                 all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
%                 all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');
% 
%                 pyr_color = [1 .7 .7];
%                 nw_color = [.7 .7 1];
%                 ww_color = [.7 1 1];
% 
%                 spikes = loadSpikes;
%                 UID = spikes.UID;
% 
%                 acg = cell_metrics.acg.log10;
%                 acg_time = cell_metrics.general.acgs.log10;
%                 acg_time_offset = acg_time(minPeakTime:end);
%                 offset = length(acg_time) - length(acg_time_offset);
% 
%                 acg_smoothed = smooth(acg,10);
%                 acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
%                 acg_smoothed_norm = acg_smoothed./sum(acg_smoothed);
%                 acg_smoothed_offset = acg_smoothed_norm(minPeakTime:end,:);
% 
%                 acgPeak = [];
%                 acgPeak2 = [];
%                 acgPeak_sample = [];
%                 acgPeak_sample2 = [];
% 
%                 for i = 1:length(UID)
%                     [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
%                     acgPeak(i) = acg_time_offset(acgPeak_sample(i));
% 
%                     [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
%                     acgPeak2(i) = acg_time(acgPeak_sample2(i));
%                 end
% 
%                 acg_time_samples = acg_time;
%                 acg_time = log10(cell_metrics.general.acgs.log10);
% 
%                 figure('position',[200 115 1300 800])
%                 % set(gcf,'Position',get(0,'screensize'));
%                 subplot(2,2,[1 2])
%                 hold on;
%                 plotFill(acg_time,acg_smoothed_norm(:,all_pyr),'Color',pyr_color);
%                 plotFill(acg_time,acg_smoothed_norm(:,all_nw),'Color',nw_color);
%                 plotFill(acg_time,acg_smoothed_norm(:,all_ww),'Color',ww_color);
% 
%                 set(gca,'XTick',[(-2) (-1) 0 1])
%                 XTick = [-2 -1 0 1];
%                 XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
%                 set(gca,'XTickLabel',XTickLabels);
%                 ylabel('logACG (prob)'); xlabel('Time(s)');
%                 axis tight;
% 
%                 subplot(2,2,3)
%                 hold on;
%                 histogram(acgPeak_sample2(all_pyr),'FaceColor',pyr_color);
%                 histogram(acgPeak_sample2(all_nw),'FaceColor',nw_color);
%                 histogram(acgPeak_sample2(all_ww),'FaceColor',ww_color);
%                 axis tight; ylabel('Count'); xlabel('bin number');xlim([0 60])
% 
%                 subplot(2,2,4)
%                 hold on;
%                 histogram(acgPeak_sample(all_pyr)+offset,'FaceColor',pyr_color);
%                 histogram(acgPeak_sample(all_nw)+offset,'FaceColor',nw_color);
%                 histogram(acgPeak_sample(all_ww)+offset,'FaceColor',ww_color);
%                 axis tight; ylabel('Count'); xlabel('bin number'); xlim([0 60])
% 
%                 saveas(gca,['BaselinevsDrug\ACGPeak_Drug.png']); 
% 
%                 acgPeak = [];
% 
%                 acgPeak.acg_smoothed = acg_smoothed;
%                 acgPeak.acg_smoothed_norm = acg_smoothed_norm;
%                 acgPeak.acgPeak_sample = acgPeak_sample+offset;
%                 acgPeak.acg_time = acg_time;
%                 acgPeak.acg_time_samples = acg_time_samples';
% 
%                 acgPeak.acgPeak_sample2 = acgPeak_sample2;
%                 acgPeak.acgPeak_sample = acgPeak_sample;
%                 acgPeak.offset = offset;
% 
%                 save([session.general.name,'.ACGPeak_Drug.cellinfo.mat'],'acgPeak');
%             else
%                 
%                 acgPeak = [];
%                 
%                 acgPeak.acg_smoothed = [];
%                 acgPeak.acg_smoothed_norm = [];
%                 acgPeak.acgPeak_sample = [];
%                 acgPeak.acg_time = [];
%                 acgPeak.acg_time_samples = [];
% 
%                 acgPeak.acgPeak_sample2 = [];
%                 acgPeak.acgPeak_sample = [];
%                 acgPeak.offset = [];
%                 
%                 save([session.general.name,'.ACGPeak_Drug.cellinfo.mat'],'acgPeak');
%             end
%         end
        
        
        % ===============================
        % 7. AVERAGE CCG
        % ==============================
        
%         if isempty(dir('*.averageCCG_Baseline.cellinfo.mat')) | isempty(dir('*.averageCCG_Drug.cellinfo.mat')) | forceACG
%             
%             winSizePlot = [-.3 .3];
%             
%             % Baseline
%             
%             averageCCG_Baseline = getAverageCCG('force',true,'includeIntervals',ts_Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%             averageCCG = averageCCG_Baseline;
%             save([session.general.name,'.averageCCG_Baseline.cellinfo.mat'],'averageCCG');
% 
%             t_ccg = averageCCG_Baseline.timestamps;
%             allCcg = averageCCG_Baseline.allCcg;
%             indCell = [1:size(allCcg,2)];
% 
%             figure('position',[200 115 1300 800])
%             for kk = 1:size(spikes.UID,2)
%                 % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                 subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                 cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                 imagesc(t_ccg,1:max(indCell)-1,cc)
%                 set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                 hold on
%                 zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                 zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                 plot(t_ccg, zmean,'k','LineWidth',2);
%                 xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                 title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                 if kk == 1
%                     ylabel('Cell');
%                 elseif kk == size(spikes.UID,2)
%                     xlabel('Time (s)');
%                 else
%                     set(gca,'YTick',[],'XTick',[]);
%                 end
%             end
%             saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Baseline.png']);
%             
%             figure('position',[200 115 1300 800])
%             imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Baseline.ZmeanCCG,1)],...
%                 averageCCG_Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
%             set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%             title('Grand CCG average','FontWeight','normal','FontSize',10);
%             
%             saveas(gca,['BaselineVsDrug\grandCCGAverage_Baseline.png']);
% 
%             % Drug 
%             
%             if ~isempty(ts_Drug)
%                 averageCCG_Drug = getAverageCCG('force',true,'includeIntervals',ts_Drug,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Drug;
%                 save([session.general.name,'.averageCCG_Drug.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Drug.timestamps;
%                 allCcg = averageCCG_Drug.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Drug.png']);
%                 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Drug.ZmeanCCG,1)],...
%                     averageCCG_Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
%                 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Drug.png']);
%             
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Drug.cellinfo.mat'],'averageCCG');
%             end
%             
%             % MazeBaseline
%             if ~isempty(ts_MazeBaseline)
%                 averageCCG_MazeBaseline = getAverageCCG('force',true,'includeIntervals',ts_MazeBaseline,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_MazeBaseline;
%                 save([session.general.name,'.averageCCG_MazeBaseline.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_MazeBaseline.timestamps;
%                 allCcg = averageCCG_MazeBaseline.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_MazeBaseline.png']);
% 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_MazeBaseline.ZmeanCCG,1)],...
%                     averageCCG_MazeBaseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
% 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_MazeBaseline.png']);
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_MazeBaseline.cellinfo.mat'],'averageCCG');
%             end
%             
%             % MazeDrug 
%             
%             if ~isempty(ts_MazeDrug)
%                 averageCCG_MazeDrug = getAverageCCG('force',true,'includeIntervals',ts_MazeDrug,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_MazeDrug;
%                 save([session.general.name,'.averageCCG_MazeDrug.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_MazeDrug.timestamps;
%                 allCcg = averageCCG_MazeDrug.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_MazeDrug.png']);
%                 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_MazeDrug.ZmeanCCG,1)],...
%                     averageCCG_MazeDrug.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
%                 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_MazeDrug.png']);
%             
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_MazeDrug.cellinfo.mat'],'averageCCG');
%             end
%             
%             % Maze1Baseline
%             if ~isempty(ts_Maze1Baseline)
%                 averageCCG_Maze1Baseline = getAverageCCG('force',true,'includeIntervals',ts_Maze1Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze1Baseline;
%                 save([session.general.name,'.averageCCG_Maze1Baseline.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze1Baseline.timestamps;
%                 allCcg = averageCCG_Maze1Baseline.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze1Baseline.png']);
% 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze1Baseline.ZmeanCCG,1)],...
%                     averageCCG_Maze1Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
% 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze1Baseline.png']);
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze1Baseline.cellinfo.mat'],'averageCCG');
%             end
%             
%             % MazeDrug 
%             
%             if ~isempty(ts_Maze1Drug)
%                 averageCCG_Maze1Drug = getAverageCCG('force',true,'includeIntervals',ts_Maze1Drug,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze1Drug;
%                 save([session.general.name,'.averageCCG_Maze1Drug.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze1Drug.timestamps;
%                 allCcg = averageCCG_Maze1Drug.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze1Drug.png']);
%                 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_MazeDrug.ZmeanCCG,1)],...
%                     averageCCG_MazeDrug.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
%                 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze1Drug.png']);
%             
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze1Drug.cellinfo.mat'],'averageCCG');
%             end
%             
%             % Maze2Baseline
%             if ~isempty(ts_Maze2Baseline)
%                 averageCCG_Maze2Baseline = getAverageCCG('force',true,'includeIntervals',ts_Maze2Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze2Baseline;
%                 save([session.general.name,'.averageCCG_Maze2Baseline.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze2Baseline.timestamps;
%                 allCcg = averageCCG_Maze2Baseline.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze2Baseline.png']);
% 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze2Baseline.ZmeanCCG,1)],...
%                     averageCCG_Maze2Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
% 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze2Baseline.png']);
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze2Baseline.cellinfo.mat'],'averageCCG');
%             end
%             
%             % Maze2Drug 
%             
%             if ~isempty(ts_Maze2Drug)
%                 averageCCG_Maze2Drug = getAverageCCG('force',true,'includeIntervals',ts_Maze2Drug,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze2Drug;
%                 save([session.general.name,'.averageCCG_Maze2Drug.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze2Drug.timestamps;
%                 allCcg = averageCCG_Maze2Drug.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze2Drug.png']);
%                 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze2Drug.ZmeanCCG,1)],...
%                     averageCCG_Maze2Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
%                 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze2Drug.png']);
%             
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze2Drug.cellinfo.mat'],'averageCCG');
%             end
%             
%             % Maze3Baseline
%             if ~isempty(ts_Maze3Baseline)
%                 averageCCG_Maze3Baseline = getAverageCCG('force',true,'includeIntervals',ts_Maze3Baseline,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze3Baseline;
%                 save([session.general.name,'.averageCCG_Maze3Baseline.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze3Baseline.timestamps;
%                 allCcg = averageCCG_Maze3Baseline.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze3Baseline.png']);
% 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze3Baseline.ZmeanCCG,1)],...
%                     averageCCG_Maze3Baseline.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
% 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze3Baseline.png']);
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze3Baseline.cellinfo.mat'],'averageCCG');
%             end
%             
%             % Maze3Drug 
%             
%             if ~isempty(ts_Maze3Drug)
%                 averageCCG_Maze3Drug = getAverageCCG('force',true,'includeIntervals',ts_Maze3Drug,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%                 averageCCG = averageCCG_Maze3Drug;
%                 save([session.general.name,'.averageCCG_Maze3Drug.cellinfo.mat'],'averageCCG');
% 
%                 t_ccg = averageCCG_Maze3Drug.timestamps;
%                 allCcg = averageCCG_Maze3Drug.allCcg;
%                 indCell = [1:size(allCcg,2)];
% 
%                 figure('position',[200 115 1300 800])
%                 for kk = 1:size(spikes.UID,2)
%                     % fprintf(' **CCG from unit %3.i/ %3.i \n',kk, size(spikes.UID,2)); %\n
%                     subplot(7,ceil(size(spikes.UID,2)/7),kk);
%                     cc = zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2); % get crosscorr
%                     imagesc(t_ccg,1:max(indCell)-1,cc)
%                     set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
%                     hold on
%                     zmean = mean(zscore(squeeze(allCcg(:,kk,indCell(indCell~=kk)))',[],2));
%                     zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
%                     plot(t_ccg, zmean,'k','LineWidth',2);
%                     xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 max(indCell)-1]);
%                     title(num2str(kk),'FontWeight','normal','FontSize',10);
% 
%                     if kk == 1
%                         ylabel('Cell');
%                     elseif kk == size(spikes.UID,2)
%                         xlabel('Time (s)');
%                     else
%                         set(gca,'YTick',[],'XTick',[]);
%                     end
%                 end
%                 saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Maze3Drug.png']);
%                 
%                 figure('position',[200 115 1300 800])
%                 imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Maze3Drug.ZmeanCCG,1)],...
%                     averageCCG_Maze3Drug.ZmeanCCG); caxis([-3 3]); colormap(jet);
%                 set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%                 title('Grand CCG average','FontWeight','normal','FontSize',10);
%                 
%                 saveas(gca,['BaselineVsDrug\grandCCGAverage_Maze3Drug.png']);
%             
%             else
%                 averageCCG = [];
%                 save([session.general.name,'.averageCCG_Maze3Drug.cellinfo.mat'],'averageCCG');
%             end
%             
%         end
        
        
        %% ====================================
        % 8. SPEED CORR
        %% ====================================
%         if isempty(dir('*.speedCorrs_Baseline.cellinfo.mat')) | isempty(dir('*.speedCorrs_Drug.cellinfo.mat')) | force
%             
%             % Baseline
%             try
%                 speedCorr_Baseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                 speedCorrs = speedCorr_Baseline;
%                 save([session.general.name,'.speedCorrs_Baseline.cellinfo.mat'],'speedCorrs');
% 
%             end
%             
%             % Drug
%             if ~isempty(ts_Drug)
%                 try
%                     speedCorr_Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Drug;
%                     save([session.general.name,'.speedCorrs_Drug.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Drug.cellinfo.mat'],'speedCorrs');
%             end
%             
%             % MazeBaseline
%             if ~isempty(ts_MazeBaseline)
%                 try
%                     speedCorr_MazeBaseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_MazeBaseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_MazeBaseline;
%                     save([session.general.name,'.speedCorrs_MazeBaseline.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_MazeBaseline.cellinfo.mat'],'speedCorrs');
%             end
% 
%             % MazeDrug
%             if ~isempty(ts_MazeDrug)
%                 try
%                     speedCorr_MazeDrug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_MazeDrug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_MazeDrug;
%                     save([session.general.name,'.speedCorrs_MazeDrug.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_MazeDrug.cellinfo.mat'],'speedCorrs');
%             end
%             
%             % Maze1Baseline
%             if ~isempty(ts_Maze1Baseline)
%                 try
%                     speedCorr_MazeBaseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze1Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze1Baseline;
%                     save([session.general.name,'.speedCorrs_Maze1Baseline.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze1Baseline.cellinfo.mat'],'speedCorrs');
%             end
% 
%             % Maze1Drug
%             if ~isempty(ts_Maze1Drug)
%                 try
%                     speedCorr_Maze1Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze1Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze1Drug;
%                     save([session.general.name,'.speedCorrs_Maze1Drug.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze1Drug.cellinfo.mat'],'speedCorrs');
%             end
%             
%             
%             % Maze2Baseline
%             if ~isempty(ts_Maze2Baseline)
%                 try
%                     speedCorr_Maze2Baseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze2Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze2Baseline;
%                     save([session.general.name,'.speedCorrs_Maze2Baseline.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze2Baseline.cellinfo.mat'],'speedCorrs');
%             end
% 
%             % Maze2Drug
%             if ~isempty(ts_Maze2Drug)
%                 try
%                     speedCorr_Maze2Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze2Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze2Drug;
%                     save([session.general.name,'.speedCorrs_Maze2Drug.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze2Drug.cellinfo.mat'],'speedCorrs');
%             end
%             
%             % Maze3Baseline
%             if ~isempty(ts_Maze3Baseline)
%                 try
%                     speedCorr_Maze3Baseline = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze3Baseline,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze3Baseline;
%                     save([session.general.name,'.speedCorrs_Maze3Baseline.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze3Baseline.cellinfo.mat'],'speedCorrs');
%             end
% 
%             % MazeDrug
%             if ~isempty(ts_Maze3Drug)
%                 try
%                     speedCorr_Maze3Drug = getSpeedCorr('numQuantiles',20,'restrictIntervals',ts_Maze3Drug,'force',true,'trials',false,'saveMat',false,'plt',false,'mkplt',false);
% 
%                     speedCorrs = speedCorr_Maze3Drug;
%                     save([session.general.name,'.speedCorrs_Maze3Drug.cellinfo.mat'],'speedCorrs');
%                 end
%             else
%                 speedCorrs = [];
%                 save([session.general.name,'.speedCorrs_Maze3Drug.cellinfo.mat'],'speedCorrs');
%             end
%             
%         end
%         
        
        %% =====================================
        % 9. SUMMARY
        %% ====================================
        
        % Baseline
%         cd('BaselinevsDrug');
%         if ~exist('Summary_Baseline.png','file')
%             cd ..
%             plotSummary_Baseline('excludePlot',{'spatialModulation'});
%         else
%             cd ..
%         end
%         % Drug
%         if ~isempty(ts_Drug)
%             cd('BaselinevsDrug')
%             if ~exist('Summary_Drug.png','file')
%                 cd ..
%                 plotSummary_Drug('excludePlot',{'spatialModulation'});
%             else 
%                 cd ..
%             end  
%         end
    end
%     close all;
    clc;
end