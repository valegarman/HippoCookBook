%% BatchScript_analysis_HM

clear; close all
targetProject= 'HMProject';
cd('J:\data');
database_path = 'J:\data';
HCB_directory = what('HMProject'); 

% sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
sessionsTable = readtable(['C:\Users\Jorge\Documents\GitHub\HMProject',filesep,'indexedSessions_HMProject.csv']);

theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

force = false;
forceACG = true;

for ii = 1:length(sessionsTable.SessionName)

    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        close all;
        % Load session info
        session = loadSession();
        spikes = loadSpikes();
        
%         cell_metrics = CellExplorer('basepath',pwd);
        
        getAverageCCG('force',true);

        ts_Pre_Spontaneous = [];
        ts_Pre_EyesClosed = [];
        ts_Pre_EyesOpen = [];
        ts_Trial = [];
        ts_interTrial = [];
        ts_Post = [];
        ts_PostEyesClosed = [];
        ts_PostEyesOpen = [];
        ts_Onomatopeyas = [];
        
        count = 1;
        for jj = 1:length(session.epochs)
            session.epochs{jj}.behavioralParadigm
            if strcmpi(session.epochs{jj}.behavioralParadigm , 'PreSleep Spontaneous')
                
                ts_Pre_Spontaneous  = [ts_Pre_Spontaneous ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm , 'PreSleep EyesClosed') 
                
                ts_Pre_EyesClosed = [ts_Pre_EyesClosed ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'PreSleep EyesOpen')
                
                ts_Pre_EyesOpen = [ts_Pre_EyesOpen ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];

            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial1') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial2') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial3') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial4') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial5')
                
                ts_Trial = [ts_Trial ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
            
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial1 Onomatopeyas') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial2 Onomatopeyas') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Trial3 Onomatopeyas') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Onomatopeyas')
    
                ts_Onomatopeyas= [ts_Onomatopeyas ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'Intertrial1') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Intertrial2') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Intertrial3') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Intertrial4') || strcmpi(session.epochs{jj}.behavioralParadigm, 'Intertrial5')
                
                ts_interTrial = [ts_interTrial ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'PostSleep')
                
                ts_Post = [ts_Post ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'PostSleep EyesOpen') || strcmpi(session.epochs{jj}.behavioralParadigm, 'PostSleep OpenEyes')
                
                ts_PostEyesOpen = [ts_PostEyesOpen ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];    
                
            elseif strcmpi(session.epochs{jj}.behavioralParadigm, 'PostSleep EyesClosed') || strcmpi(session.epochs{jj}.behavioralParadigm, 'PostSleep ClosedEyes')
                
                ts_PostEyesClosed = [ts_PostEyesClosed ; session.epochs{jj}.startTime session.epochs{jj}.stopTime];      
                
            end
        end
        
        %% ========================================================
        % 1. AVERAGE CCG 
        % =========================================================
        
        if isempty(dir('*.averageCCG_PreSpontaneous.cellinfo.mat')) || forceACG
            
            winSizePlot = [-.3 .3];
            
            % Pre Spontaneous
            if ~isempty(ts_Pre_Spontaneous)
                averageCCG_Pre_Spontaneous = getAverageCCG('force',true,'includeIntervals',ts_Pre_Spontaneous,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Pre_Spontaneous;
                save([session.general.name,'.averageCCG_PreSpontaneous.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Pre_Spontaneous.timestamps;
                allCcg = averageCCG_Pre_Spontaneous.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_PreSpontaneous.png']);

                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Pre_Spontaneous.ZmeanCCG,1)],...
                    averageCCG_Pre_Spontaneous.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);

                saveas(gca,['BaselineVsDrug\grandCCGAverage_PreSpontaneous.png']);
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_PreSpontaneous.cellinfo.mat'],'averageCCG');
            end
            
            
            % Pre Eyes Closed
            if ~isempty(ts_Pre_EyesClosed)
                averageCCG_Pre_EyesClosed = getAverageCCG('force',true,'includeIntervals',ts_Pre_EyesClosed,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Pre_EyesClosed;
                save([session.general.name,'.averageCCG_PreEyesClosed.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Pre_EyesClosed.timestamps;
                allCcg = averageCCG_Pre_EyesClosed.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_PreEyesClosed.png']);

                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Pre_EyesClosed.ZmeanCCG,1)],...
                    averageCCG_Pre_EyesClosed.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);

                saveas(gca,['BaselineVsDrug\grandCCGAverage_PreEyesClosed.png']);
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_PreEyesClosed.cellinfo.mat'],'averageCCG');
            end

            % Pre Eyes Open
            
            if ~isempty(ts_Pre_EyesOpen)
                averageCCG_Pre_EyesOpen = getAverageCCG('force',true,'includeIntervals',ts_Pre_EyesOpen,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Pre_EyesOpen;
                save([session.general.name,'.averageCCG_PreEyesOpen.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Pre_EyesOpen.timestamps;
                allCcg = averageCCG_Pre_EyesOpen.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_PreEyesOpen.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Pre_EyesOpen.ZmeanCCG,1)],...
                    averageCCG_Pre_EyesOpen.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_PreEyesOpen.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_PreEyesOpen.cellinfo.mat'],'averageCCG');
            end
            
            % Trial
            
            if ~isempty(ts_Trial)
                averageCCG_Trial = getAverageCCG('force',true,'includeIntervals',ts_Trial,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Trial;
                save([session.general.name,'.averageCCG_Trial.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Trial.timestamps;
                allCcg = averageCCG_Trial.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Trial.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Trial.ZmeanCCG,1)],...
                    averageCCG_Trial.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_Trial.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_Trial.cellinfo.mat'],'averageCCG');
            end
            
            % interTrial
            
            if ~isempty(ts_interTrial)
                averageCCG_interTrial = getAverageCCG('force',true,'includeIntervals',ts_interTrial,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_interTrial;
                save([session.general.name,'.averageCCG_interTrial.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_interTrial.timestamps;
                allCcg = averageCCG_interTrial.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_interTrial.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_interTrial.ZmeanCCG,1)],...
                    averageCCG_interTrial.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_interTrial.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_interTrial.cellinfo.mat'],'averageCCG');
            end
            
            
            % Post
            
            if ~isempty(ts_Post)
                averageCCG_Post = getAverageCCG('force',true,'includeIntervals',ts_Post,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_Post;
                save([session.general.name,'.averageCCG_Post.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_Post.timestamps;
                allCcg = averageCCG_Post.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_Post.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_Post.ZmeanCCG,1)],...
                    averageCCG_Post.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_Post.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_Post.cellinfo.mat'],'averageCCG');
            end
            
            % PostEyesClosed
            
            if ~isempty(ts_PostEyesClosed)
                averageCCG_PostEyesClosed = getAverageCCG('force',true,'includeIntervals',ts_PostEyesClosed,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_PostEyesClosed;
                save([session.general.name,'.averageCCG_PostEyesClosed.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_PostEyesClosed.timestamps;
                allCcg = averageCCG_PostEyesClosed.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_PostEyesClosed.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_PostEyesClosed.ZmeanCCG,1)],...
                    averageCCG_PostEyesClosed.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_PostEyesClosed.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_PostEyesClosed.cellinfo.mat'],'averageCCG');
            end
            
            % Post Eyes Open
            if ~isempty(ts_PostEyesOpen)
                averageCCG_PostEyesOpen = getAverageCCG('force',true,'includeIntervals',ts_PostEyesOpen,'savemat',false,'plotOpt',false,'saveFig',false);

                averageCCG = averageCCG_PostEyesOpen;
                save([session.general.name,'.averageCCG_PostEyesOpen.cellinfo.mat'],'averageCCG');

                t_ccg = averageCCG_PostEyesOpen.timestamps;
                allCcg = averageCCG_PostEyesOpen.allCcg;
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
                saveas(gca,['BaselineVsDrug\allCellsAverageCCG_PostEyesOpen.png']);
                
                figure('position',[200 115 1300 800])
                imagesc([t_ccg(1) t_ccg(end)],[1 size(averageCCG_PostEyesOpen.ZmeanCCG,1)],...
                    averageCCG_PostEyesOpen.ZmeanCCG); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([winSizePlot(1) winSizePlot(2)]);
                title('Grand CCG average','FontWeight','normal','FontSize',10);
                
                saveas(gca,['BaselineVsDrug\grandCCGAverage_PostEyesOpen.png']);
            
            else
                averageCCG = [];
                save([session.general.name,'.averageCCG_PostEyesOpen.cellinfo.mat'],'averageCCG');
            end
            
        end
        
        % ===============================
        % 5. CELL METRICS        
        % ===============================
        
        if isempty(dir('*.cell_metrics_PreSpontaneous.cellinfo.mat')) | force
            
            if isempty(ts_Pre_EyesClosed)
                
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_PreSpontaneous.cellinfo.mat'],'cell_metrics');
            else
                
                cell_metrics_Spontaneous = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Pre_Spontaneous,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_Spontaneous;
                save([session.general.name,'.cell_metrics_PreSpontaneous.cellinfo.mat'],'cell_metrics');
            end
           
        end
        
        if isempty(dir('*.cell_metrics_PreEyesClosed.cellinfo.mat')) | force
            
            if isempty(ts_Pre_EyesClosed)
                
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_PreEyesClosed.cellinfo.mat'],'cell_metrics');
            else
                
                cell_metrics_PreEyesClosed = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Pre_EyesClosed,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_PreEyesClosed;
                save([session.general.name,'.cell_metrics_PreEyesClosed.cellinfo.mat'],'cell_metrics');
            end
           
        end

        
        if isempty(dir('*.cell_metrics_PreEyesOpen.cellinfo.mat')) | force
            
            if isempty(ts_Pre_EyesOpen)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_PreEyesOpen.cellinfo.mat'],'cell_metrics');
            else
            
                cell_metrics_PreEyesOpen = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Pre_EyesOpen,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_PreEyesOpen;
                save([session.general.name,'.cell_metrics_PreEyesOpen.cellinfo.mat'],'cell_metrics');
            end
        end
        
        if isempty(dir('*.cell_metrics_Trial.cellinfo.mat')) | force  
            if isempty(ts_Trial)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_Trial.cellinfo.mat'],'cell_metrics');
            else
            
                cell_metrics_Trial = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Trial,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_Trial;
                save([session.general.name,'.cell_metrics_Trial.cellinfo.mat'],'cell_metrics');
            end 
        end
        
        if isempty(dir('*.cell_metrics_Onomatopeyas.cellinfo.mat')) | force  
            if isempty(ts_Onomatopeyas)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_Onomatopeyas.cellinfo.mat'],'cell_metrics');
            else
            
                cell_metrics_Onomatopeyas = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Onomatopeyas,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_Onomatopeyas;
                save([session.general.name,'.cell_metrics_Onomatopeyas.cellinfo.mat'],'cell_metrics');
            end 
        end
        
        if isempty(dir('*.cell_metrics_interTrial.cellinfo.mat')) | force  
            if isempty(ts_interTrial)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_interTrial.cellinfo.mat'],'cell_metrics');
            else
                cell_metrics_interTrial = ProcessCellMetrics('session', session,'restrictToIntervals',ts_interTrial,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_interTrial;
                save([session.general.name,'.cell_metrics_interTrial.cellinfo.mat'],'cell_metrics');
            end
           
        end
        
        if isempty(dir('*.cell_metrics_Post.cellinfo.mat')) | force
            if isempty(ts_Post)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_Post.cellinfo.mat'],'cell_metrics');
            else
            
                cell_metrics_Post = ProcessCellMetrics('session', session,'restrictToIntervals',ts_Post,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
                cell_metrics = cell_metrics_Post;
                save([session.general.name,'.cell_metrics_Post.cellinfo.mat'],'cell_metrics');
            end
        end
        
        if isempty(dir('*.cell_metrics_PostEyesOpen.cellinfo.mat')) | force
            if isempty(ts_PostEyesOpen)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_PostEyesOpen.cellinfo.mat'],'cell_metrics');
            else
                cell_metrics_PostEyesOpen = ProcessCellMetrics('session', session,'restrictToIntervals',ts_PostEyesOpen,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);
                cell_metrics = cell_metrics_PostEyesOpen;
                save([session.general.name,'.cell_metrics_PostEyesOpen.cellinfo.mat'],'cell_metrics');
            end
        end
        
        if isempty(dir('*.cell_metrics_PostEyesClosed.cellinfo.mat')) | force
            if isempty(ts_PostEyesClosed)
                cell_metrics = [];
                save([session.general.name,'.cell_metrics_PostEyesClosed.cellinfo.mat'],'cell_metrics');
                
            else
                cell_metrics_PostEyesClosed = ProcessCellMetrics('session', session,'restrictToIntervals',ts_PostEyesClosed,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'},'forceReload',true,'saveMat',false);

                cell_metrics = cell_metrics_PostEyesClosed;
                save([session.general.name,'.cell_metrics_PostEyesClosed.cellinfo.mat'],'cell_metrics');
            end
           
        end
        
        % =================================
        % 6. PHASE MODULATION
        % ================================
        
        if isempty(dir('*theta_6-12_PreSpontaneous.PhaseLockingData.cellinfo.mat')) || force
            try
                targetFile = dir('*theta_6-12.PhaseLockingData.cellinfo.mat');
                load(targetFile.name);
                thetaChannel = thetaMod.lfpChannel;

            catch
                disp('No thetaChannel detected!');
            end

            rippleChannel = [];
            SWChannel = [];
            mkdir('BaselinevsDrug')

            % PreSpontaneous
            if isempty(ts_Pre_Spontaneous)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PreSpontaneous.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '_PreSpontaneous.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PreSpontaneous.PhaseLockingData.cellinfo.mat'],'hgammaMod');

            else

                phaseMod_PreSpontaneous = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Pre_Spontaneous,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_PreSpontaneous.theta;
                lgammaMod = phaseMod_PreSpontaneous.lgamma;
                hgammaMod = phaseMod_PreSpontaneous.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PreSpontaneous.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PreSpontaneous_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_PreEyesSpontaneous.PhaseLockingData.cellinfo.mat'],'lgammaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PreSpontaneous_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PreSpontaneous.PhaseLockingData.cellinfo.mat'],'hgammaMod');   
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PreSpontaneous_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                 end
            end


            % PreEyesClosed
            if isempty(ts_Pre_EyesClosed)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'hgammaMod');

            else

                phaseMod_PreEyesClosed = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Pre_EyesClosed,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_PreEyesClosed.theta;
                lgammaMod = phaseMod_PreEyesClosed.lgamma;
                hgammaMod = phaseMod_PreEyesClosed.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PreEyesClosed_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'lgammaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PreEyesClosed_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PreEyesClosed.PhaseLockingData.cellinfo.mat'],'hgammaMod');   
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PreEyesClosed_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                 end
            end

            % PreEyesOpen
            if isempty(ts_Pre_EyesOpen)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];
                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'hgammaMod');

            else
                phaseMod_PreEyesOpen = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Pre_EyesOpen,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_PreEyesOpen.theta;
                lgammaMod = phaseMod_PreEyesOpen.lgamma;
                hgammaMod = phaseMod_PreEyesOpen.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PreEyesOpen_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'lgammaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PreEyesOpen_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PreEyesOpen.PhaseLockingData.cellinfo.mat'],'hgammaMod');   

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PreEyesOpen_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end
            end

           % Trial
           if isempty(ts_Trial)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Trial.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Trial.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Trial.PhaseLockingData.cellinfo.mat'],'hgammaMod');
           else

                phaseMod_Trial = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Trial,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_Trial.theta;
                lgammaMod = phaseMod_Trial.lgamma;
                hgammaMod = phaseMod_Trial.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_Trial.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Trial_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_Trial.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Trial_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_Trial.PhaseLockingData.cellinfo.mat'],'hgammaMod');  
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Trial_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end
           end


           % Onomatopeyas
           if isempty(ts_Onomatopeyas)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'hgammaMod');
           else

                phaseMod_Onomatopeyas = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Onomatopeyas,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_Onomatopeyas.theta;
                lgammaMod = phaseMod_Onomatopeyas.lgamma;
                hgammaMod = phaseMod_Onomatopeyas.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Onomatopeyas_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Onomatopeyas_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                 save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_Onomatopeyas.PhaseLockingData.cellinfo.mat'],'hgammaMod');  
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Onomatopeyas_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end
           end

            % InterTrial
            if isempty(ts_interTrial)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_interTrial.PhaseLockingData.cellinfo.mat'],'thetaMod');

            else

                phaseMod_interTrial = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_interTrial,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_interTrial.theta;
                lgammaMod = phaseMod_interTrial.lgamma;
                hgammaMod = phaseMod_interTrial.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_interTrial.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_interTrial.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_interTrial.PhaseLockingData.cellinfo.mat'],'hgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_interTrial_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_interTrial.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_interTrial_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_interTrial.PhaseLockingData.cellinfo.mat'],'hgammaMod');  
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_interTrial_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end
            end

            % Post
            if isempty(ts_Post)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_Post.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_Post.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_Post.PhaseLockingData.cellinfo.mat'],'hgammaMod');
            else     
                phaseMod_Post = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_Post,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_Post.theta;
                lgammaMod = phaseMod_Post.lgamma;
                hgammaMod = phaseMod_Post.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_Post.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_Post_PhaseModulation.png']);
                   catch
                       disp('Not possible to run theta modulation Baseline vs Drug...');
                    end

                    save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                            '_Post.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_Post_PhaseModulation.png'])
                    catch
                       disp('Not possible to run lgamma modulation Baseline vs Drug...');
                    end

                    save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                            '_Post.PhaseLockingData.cellinfo.mat'],'hgammaMod');  

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
                                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_Post_PhaseModulation.png']);
                   catch
                       disp('Not possible to run hgamma modulation Baseline vs Drug...');
                    end    
                end

            % PostEyesOpen
            if isempty(ts_PostEyesOpen)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'hgammaMod');
            else     
                phaseMod_PostEyesOpen = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_PostEyesOpen,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_PostEyesOpen.theta;
                lgammaMod = phaseMod_PostEyesOpen.lgamma;
                hgammaMod = phaseMod_PostEyesOpen.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PostEyesOpen_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PostEyesOpen_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PostEyesOpen.PhaseLockingData.cellinfo.mat'],'hgammaMod');  

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PostEyesOpen_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end    
            end

            % PostEyesClosed

            if isempty(ts_PostEyesClosed)
                thetaMod = [];
                lgammaMod = [];
                hgammaMod = [];

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'thetaMod');
                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                    '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'lgammaMod');
                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                    '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'hgammaMod');
            else     
                phaseMod_PostEyesClosed = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrictIntervals',ts_PostEyesClosed,'plotting',false,'saveMat',false);

                thetaMod = phaseMod_PostEyesClosed.theta;
                lgammaMod = phaseMod_PostEyesClosed.lgamma;
                hgammaMod = phaseMod_PostEyesClosed.hgamma;

                save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                        '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'thetaMod');

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PostEyesClosed_PhaseModulation.png']);
               catch
                   disp('Not possible to run theta modulation Baseline vs Drug...');
                end

                save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                        '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'lgammaMod');
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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PostEyesClosed_PhaseModulation.png'])
                catch
                   disp('Not possible to run lgamma modulation Baseline vs Drug...');
                end

                save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                        '_PostEyesClosed.PhaseLockingData.cellinfo.mat'],'hgammaMod');  

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
                            ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
                        elseif i == size(spikes.UID,2)
                            set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                            xlabel('phase (rad)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gca,['BaselinevsDrug\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PostEyesClosed_PhaseModulation.png']);
               catch
                   disp('Not possible to run hgamma modulation Baseline vs Drug...');
                end    
            end
        end
        
    end
    
    close all;
    clc;
    
end