%% BatchScript_analysis_modulationPerSubSession

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
forceReload = false;
plt = true;
ripple_passband = [120 200];
SW_passband = [2 10];
theta_passband = [6 12];
lgamma_passband = [20 60];
hgamma_passband = [60 100];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        if isempty(dir(['*theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'.PhaseLockingData.Subsession.cellinfo.mat'])) | forceReload
            targetFile = dir('*thetaEpochs.states.mat');
            try
                load(targetFile.name);
            catch
                warning('thetaEpochs not found. Quitting analysis...');
            end
            targetFile = dir('*.ripples.events.mat');
            try
                load(targetFile.name);
            catch
                warning('ripples not found. Quitting analysis...');
            end
            
            spikes = loadSpikes;
            
            % Load session info
            session = loadSession();
            for jj = 1:length(session.epochs)
                ts = [session.epochs{jj}.startTime session.epochs{jj}.stopTime];
                
                try
                    [phaseMod] = computePhaseModulation('rippleChannel',session.analysisTags.rippleChannel,'SWChannel',session.analysisTags.SWChannel,'thetaChannel',session.analysisTags.thetaChannel,'lgammaChannel',session.analysisTags.thetaChannel,'hgammaChannel',session.analysisTags.thetaChannel,...
                                        'restrictIntervals',ts,'saveMat',false,'plotting',false);
                    
                    % Saving in subfolder
                    rippleMod.(session.epochs{jj}.name) = phaseMod.ripples;
                    SWMod.(session.epochs{jj}.name) = phaseMod.SharpWave;
                    thetaMod.(session.epochs{jj}.name) = phaseMod.theta;
                    lgammaMod.(session.epochs{jj}.name) = phaseMod.lgamma;
                    hgammaMod.(session.epochs{jj}.name) = phaseMod.hgamma;
                    thetaRunMod.(session.epochs{jj}.name) = phaseMod.thetaRunMod;
                    thetaREMMod.(session.epochs{jj}.name) = phaseMod.thetaREMMod;
                                        
                catch
                    rippleMod.(session.epochs{jj}.name) = [];
                    SWMod.(session.epochs{jj}.name) = [];
                    thetaMod.(session.epochs{jj}.name) = [];
                    lgammaMod.(session.epochs{jj}.name) = [];
                    hgammaMod.(session.epochs{jj}.name) = [];
                    thetaRunMod.(session.epochs{jj}.name) = [];
                    thetaREMMod.(session.epochs{jj}.name) = [];
                end
                
                % Saving variables
                
                    
                if plt
                    % Ripple Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.ripples.phasebins ; phaseMod.ripples.phasebins + 2*pi],[phaseMod.ripples.phasedistros(:,i) ;  phaseMod.ripples.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(session.analysisTags.rippleChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                    catch
                    end
                    
                    % SharpWave Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.SharpWave.phasebins ; phaseMod.SharpWave.phasebins + 2*pi],[phaseMod.SharpWave.phasedistros(:,i) ;  phaseMod.SharpWave.phasedistros(:,i)], 'EdgeColor','none');
                            hold on;
                            ax = axis;
                            x = 0:.001:4*pi;
                            y = cos(x);
                            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                            xlim([0 4*pi]);
                            title(num2str(i),'FontWeight','normal','FontSize',10);
                            if i == 1
                                ylabel('prob'); title(['Channel (1-index): ' num2str(session.analysisTags.SWChannel)],'FontWeight','normal','FontSize',10);
                            elseif i == size(spikes.UID,2)
                                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                xlabel('phase (rad)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                    catch
                    end
                    
                    % Theta Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.theta.phasebins ; phaseMod.theta.phasebins + 2*pi],[phaseMod.theta.phasedistros(:,i) ;  phaseMod.theta.phasedistros(:,i)], 'EdgeColor','none');
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
                    catch
                    end
                    
                    % Low gamma modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.lgamma.phasebins ; phaseMod.lgamma.phasebins + 2*pi],[phaseMod.lgamma.phasedistros(:,i) ;  phaseMod.lgamma.phasedistros(:,i)], 'EdgeColor','none');
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
                    catch
                    end
                    
                    % High Gamma Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.hgamma.phasebins ; phaseMod.hgamma.phasebins + 2*pi],[phaseMod.hgamma.phasedistros(:,i) ;  phaseMod.hgamma.phasedistros(:,i)], 'EdgeColor','none');
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
                    catch
                    end
                    % Theta RUN Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.thetaRunMod.phasebins ; phaseMod.thetaRunMod.phasebins + 2*pi],[phaseMod.thetaRunMod.phasedistros(:,i) ;  phaseMod.thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
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
                    catch
                    end
                    
                    % ThetaREM Modulation
                    try
                        figure,
                        set(gcf,'Position',get(0,'ScreenSize'))
                        for i = 1:length(spikes.UID)
                            subplot(7,ceil(size(spikes.UID,2)/7),i)
                            area([phaseMod.thetaREMMod.phasebins ; phaseMod.thetaREMMod.phasebins + 2*pi],[phaseMod.thetaREMMod.phasedistros(:,i) ;  phaseMod.thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
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
                    catch
                    end
                    
                end
            end
            
            % Ripples
            save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
                    '.PhaseLockingData.Subsession.cellinfo.mat'],'rippleMod');
            % Sharp Wave
            save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'SWMod');
            % Theta
            save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'thetaMod');
            % Low Gamma
            save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'lgammaMod');
            % High Gamma
            save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'hgammaMod');
            % ThetaRun
            save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'thetaRunMod');
            % thetaREM
            save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                '.PhaseLockingData.Subsession.cellinfo.mat'],'thetaREMMod');   
        end
    end
    close all;
end
