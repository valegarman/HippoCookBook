%% BatchScript_analysis_ripplesPerSubsession

clear; close all
targetProject= 'SubiculumProject';
cd('D:\');
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

for ii = 6:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
%         getAverageCCG('force',true);
        for jj = 1:length(session.epochs)
            
            ts = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];  
            
            % ======================
            % 1. CELL METRICS 
            % ======================
            
%             try
%                 if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
%                     file = dir([session.general.name,'.optogeneticPulses.events.mat']);
%                     load(file.name);
%                 end
%                 excludeManipulationIntervals = optoPulses.stimulationEpochs;
%             catch
%                 warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
%                 excludeManipulationIntervals = [];
%             end
% 
%             cell_metrics_epoch = ProcessCellMetrics('session', session,'restrictToIntervals',ts,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
% 
%             cell_metrics = cell_metrics_epoch;
%             save([session.general.name,'.cell_metrics_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'cell_metrics');
            
            % ========================
            % 2. CCG
            % ========================
            
%             winSizePlot = [-.3 .3];
% 
%             averageCCG_epoch = getAverageCCG('force',true,'includeIntervals',ts,'savemat',false,'plotOpt',false,'saveFig',false);
% 
%             averageCCG = averageCCG_epoch;
%             save([session.general.name,'.averageCCG_',session.epochs{jj}.behavioralParadigm,'.cellinfo.mat'],'averageCCG');
% 
%             t_ccg = averageCCG_epoch.timestamps;
%             allCcg = averageCCG_epoch.allCcg;
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
%             saveas(gca,['Epochs\allCellsAverageCCG_',session.epochs{jj}.behavioralParadigm,'.png']);
% 
%             figure('position',[200 115 1300 800])
%             imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_epoch.ZmeanCCG,1)],...
%                 averageCCG_epoch.ZmeanCCG); caxis([-3 3]); colormap(jet);
%             set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
%             title('Grand CCG average','FontWeight','normal','FontSize',10);
% 
%             saveas(gca,['Epochs\grandCCGAverage_',session.epochs{jj}.behavioralParadigm,'.png']);
        end
        
        % ===========================
        % 3. SPATIAL MODULATION
        % ===========================
        speedThresh = 1; % 5cm/s
        pixelsPerCm = 2.5;
        randomization = true;
        gridAnalysis = false;
        twoHalvesAnalysis = true;
        
        spikes = loadSpikes;
        tracking = getSessionTracking;
        behavior = getSessionBehavior;
        firingMaps = firingMapAvg_pablo(behavior,spikes,'speedThresh',speedThresh,'tint',false,'pixelsPerCm',pixelsPerCm,'saveMat',true);
        firingMaps_tint = firingMapAvg_pablo(behavior,spikes,'speedThresh',speedThresh,'tint',true,'pixelsPerCm',pixelsPerCm,'saveMat',true);
        placeFieldStats = computeFindPlaceFields('firingMaps',[],'useColorBar',false,'saveMat',true);
        spatialModulation = computeSpatialModulation('force',true,'tint',false,'gridAnalysis',gridAnalysis,'randomization',randomization,'speedThresh',speedThresh);
        spatialModulation_tint = computeSpatialModulation('force',true,'tint',true,'gridAnalysis',gridAnalysis,'randomization',randomization,'speedThresh',speedThresh);

        if twoHalvesAnalysis
            firingMaps2Halves = firingMap2Halves(behavior,spikes,'pixelsPerCm',pixelsPerCm,'speedThresh',speedThresh,'saveMat',true,'tint',false);
            firingMaps2Halves_tint = firingMap2Halves(behavior,spikes,'pixelsPerCm',pixelsPerCm,'speedThresh',speedThresh,'saveMat',true,'tint',true);
            placeFieldStats2Halves = computeFindPlaceFields2Halves('firingMaps',[],'useColorBar',false);
            spatialModulation2Halves = computeSpatialModulation2Halves('force',true,'tint',false,'speedThresh',speedThresh); % Not running gridAnalysis for two Halves
            spatialModulation2Halves_tint = computeSpatialModulation2Halves('force',true,'tint',true,'speedThresh',speedThresh); % Not running gridAnalysis for two Halves
        end
        
        plotSpatialModulation('gridAnalysis',gridAnalysis,'tint',false);
        plotSpatialModulation('gridAnalysis',gridAnalysis,'tint',true);

    end
    close all;
end