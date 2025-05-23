
function [armChoice] = getArmChoice(varargin)
% Compute T maze performance and gets timestamps
%
% USAGE
%
%   [armChoice, behavior] = getArmChoice(varargin)
%
% INPUTS
% 
% digitalIn                     digitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left Alternation,  4.Righ Alternation
%                                 5. Home Delay,        6. Is alternation forzed?
%                                 7. 0 is cue left, 1 is cue right (for
%                                       cueSide task)
%                               NOTE: After Maze reconstruction, convention
%                               changed (09/12/2022):
%                                 1. Basler             2. Left alternation
%                                 3. Right alternation  4. Home delay
%                                 5. Opto
%                               If not called, look for it.
% task                          'alternation' and 'cudeSide'
% force                         Force detection (boolean, default false)
% verbose                       default false
%
% OUTPUT
%       - armChoice.behaviour output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.visitedArm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.timestamps    Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.choice              Performance vector, 1 is right choice, 0 is
%                                   wrong. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
% armChoice.task                'alternation' and 'cudeSide'  
% armChoice.expectedArm          In 'cueSide', it specifias the right choice. 
%
%  Manuel Valero 2019
%
%
% Modified by Pablo Abad 2022 to include performance for trials with
% different homeDelay and plotting them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'digitalIn',[],@ischar);
addParameter(p,'task',[]);
addParameter(p,'force',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'leftArmTtl_channel',2,@isnumeric);
addParameter(p,'rightArmTtl_channel',3,@isnumeric);
addParameter(p,'homeDelayTtl_channel',4,@isnumeric);
addParameter(p,'delayPlot',true,@islogical);
addParameter(p,'use_manual_ttls',false,@islogical);
addParameter(p,'save_original_ttls',true,@islogical);
addParameter(p,'dlcTtl_channel',5,@isnumeric);
addParameter(p,'fiberTtl_channel',1,@isnmueric);


parse(p,varargin{:});

digitalIn = p.Results.digitalIn;
task = p.Results.task;
force = p.Results.force;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;
delayPlot = p.Results.delayPlot;
use_manual_ttls = p.Results.use_manual_ttls;
save_original_ttls = p.Results.save_original_ttls;
dlcTtl_channel = p.Results.dlcTtl_channel;
fiberTtl_channel = p.Results.fiberTtl_channel;

if ~isempty(dir('*.ArmChoice.Events.mat')) && ~force 
    disp('Arm choice already computed! Loading file.');
    file = dir('*.ArmChoice.Events.mat');
    load(file.name);
    return
end
%% Get basler TTL
disp('Loading digital In...');
if isempty(digitalIn)
    digitalIn = getDigitalIn;
    if isempty(digitalIn)
        armChoice = [];
        return
    end
end
if isempty(task)
    switch size(digitalIn.timestampsOn,2)
        case 5
            task = 'alternation';
        case 6
            task = 'alternation';
        case 7
            task = 'cueSide';
        case 10
            task = 'alternation';
        case 16
            task = 'alternation';
    end
end
if strcmpi(task,'cueSide')
    if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
        armChoice.timestamps = [digitalIn.timestampsOn{3}; digitalIn.timestampsOn{4}]; 
        % 0 is left, 1 is right
        armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{3})); ones(size(digitalIn.timestampsOn{4}))];
        [armChoice.timestamps, idx] = sort(armChoice.timestamps);
        armChoice.visitedArm = armChoice.visitedArm(idx);
        % armChoice.delay.ints = digitalIn.ints{5};
        armChoice.delay.timestamps = digitalIn.ints{5};
        armChoice.delay.dur = nanmean(armChoice.delay.timestamps(2,:) - armChoice.delay.timestamps(1,:));
        
        % reconstruct side vector
        sideVector = zeros(1,int32(max(armChoice.timestamps)));
        onPulses = digitalIn.timestampsOn{7};
        onPulses(:,2) = nan;
        offPulses = digitalIn.timestampsOff{7};
        if  offPulses(1) < onPulses(1,1)
            offPulses(1) = [];
        end
        onPulses(1:length(offPulses),2) = offPulses;
        onPulses(isnan(onPulses)) = armChoice.timestamps(end) + 10;      
      
        for ii = 1:size(onPulses,1)
            sideVector(round(onPulses(ii,1)):round(onPulses(ii,2))) = 1;
        end
        armChoice.expectedArm = sideVector(round(armChoice.timestamps))'; % real choice;
        armChoice.choice = sideVector(round(armChoice.timestamps))' == armChoice.visitedArm; % 1 is right, 0 is wrong
        armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice));

        if size(digitalIn.timestampsOn,2) >=6
            if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
                armChoice.forzed = 1;
                desc = 'forzed side';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous side';
            end
        else
            if armChoice.performance == 1 
                armChoice.forzed = 1;
                desc = 'forzed side';
            else
                armChoice.forzed = 0;
                desc = 'spontaneous side';
            end
        end

        h = figure;
        subplot(2,2,[1 2])
        hold on
        plot(armChoice.timestamps, armChoice.expectedArm,'color',[.7 .7 .7]);
        scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
            armChoice.expectedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled','MarkerFaceAlpha',.5);
        scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
            armChoice.expectedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled','MarkerFaceAlpha',.5);
        for ii = 1:size(armChoice.delay.timestamps,1)
            fill([armChoice.delay.timestamps(:,ii); flip(armChoice.delay.timestamps(:,ii))],[1 1 1.2 1.2]',...
            [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
        end
        xlabel('seconds'); ylim([-.2 1.2]);
        text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
            desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
            num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
            's'));
        set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});
        subplot(2,2,3)
        plot(cumsum(~armChoice.choice));
        set(gca,'TickDir','out'); xlabel('Trials'); ylabel('Errors');
 
        subplot(2,2,4)
        hold on
        for ii = 1:length(armChoice.choice)- 9
            lcurve(ii)=sum(armChoice.choice(ii:ii+9))/10;
        end
        lcurve_m = [];
        g = 1:10:length(armChoice.choice);
        for ii = 1:length(g)
            if g(ii) + 10 <= length(armChoice.choice)
                lcurve_m(ii) = mean(armChoice.choice(g(ii):g(ii)+10));
            else
                lcurve_m(ii) = mean(armChoice.choice(g(ii):length(armChoice.choice)));
            end
        end
        ecurve = cumsum(armChoice.choice)./(1:length(armChoice.choice))';
        area(1:length(armChoice.choice),ecurve,'FaceColor',[.8 .5 1],'EdgeColor','none');
        scatter(1:10:length(armChoice.choice),lcurve_m,100,[.9 .6 .7],'filled','MarkerFaceAlpha',.5);
        ylim([0 1]); xlabel('Trials'); ylabel('Performance'); 
        
        armChoice.lcurve_m = lcurve_m;
        armChoice.task = task;

        mkdir('Behavior');
        saveas(h,'Behavior\armChoice.png');
        if ~verbose
            close(h);
        end

        C = strsplit(pwd,'\');
        save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
    else
        warning('DigitalIn format does not match. Was Cue Side Maze performed? ');
    end
elseif strcmpi(task,'alternation')

    if ~use_manual_ttls
        % score for alternation task
        if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
            armChoice.timestamps = [digitalIn.timestampsOn{leftArmTtl_channel}; digitalIn.timestampsOn{rightArmTtl_channel}]; 
            % 0 is left, 1 is right
            armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{leftArmTtl_channel})); ones(size(digitalIn.timestampsOn{rightArmTtl_channel}))];
            armChoice.delay.timestamps = digitalIn.ints{homeDelayTtl_channel};
            if size(armChoice.visitedArm,1) < size(digitalIn.timestampsOn{homeDelayTtl_channel},1) - 10
                warning('There was problem with one of the sensors! Triying to fix it!')
                leftArmSensor = digitalIn.timestampsOn{leftArmTtl_channel};
                rightArmSensor = digitalIn.timestampsOn{rightArmTtl_channel};
                delaySensor = digitalIn.timestampsOn{homeDelayTtl_channel};
                if size(leftArmSensor,1) < size(rightArmSensor,1) -10
                    warning('Estimating left-arm sensor times...');
                    for ii = 2:length(delaySensor)-1
                        if any(leftArmSensor > delaySensor(ii) & leftArmSensor < delaySensor(ii+1)) | ...
                                any(rightArmSensor > delaySensor(ii) & rightArmSensor < delaySensor(ii+1))
                        else
                            leftArmSensor = [leftArmSensor; mean(delaySensor(ii:ii+1))];
                        end
                    end
                elseif size(rightArmSensor,1) < size(leftArmSensor,1) -10
                    warning('Estimating left-arm sensor times...');
                    for ii = 2:length(delaySensor)-1
                        if any(leftArmSensor > delaySensor(ii) & leftArmSensor < delaySensor(ii+1)) | ...
                                any(rightArmSensor > delaySensor(ii) & rightArmSensor < delaySensor(ii+1))
                        else
                            rightArmSensor = [rightArmSensor; mean(delaySensor(ii:ii+1))];
                        end
                    end
                else error('Not able to guess what is the missing sensor!!');
                end
                armChoice.timestamps = [leftArmSensor; rightArmSensor];
                armChoice.visitedArm = [zeros(size(leftArmSensor)); ones(size(rightArmSensor))];
            end
            
            [armChoice.timestamps, idx] = sort(armChoice.timestamps);
            armChoice.visitedArm = armChoice.visitedArm(idx);
            
            armChoice.delay.dur = nanmean(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1));
            armChoice.delay.durations = round(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1),0);
            armChoice.choice = [NaN; abs(diff(armChoice.visitedArm))]; % 1 is right, 0 is wrong
            armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice) - 1);
            performance = [];
            
            if ~isnan(armChoice.delay.durations(end))
                armChoice.delay.durations(end) = NaN;
            end
            durations = unique(armChoice.delay.durations);
            for ii = 1:length(durations)- 1
                performance = [performance; sum(armChoice.choice(find(armChoice.delay.durations == durations(ii)) + 1)) / length(find(armChoice.delay.durations == durations(ii)))];
            end
            armChoice.delay.performance = performance;
            armChoice.delay.uniqueDurations = durations;
            armChoice.task = task;
            armChoice.expectedArm = [NaN; ~xor(armChoice.visitedArm(2:end), armChoice.choice(2:end))];
    
            if size(digitalIn.timestampsOn,2) >=6
                if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
                    armChoice.forzed = 1;
                    desc = 'forzed alternation';
                else
                    armChoice.forzed = 0;
                    desc = 'spontaneous alternation';
                end
            else
                if armChoice.performance == 1 
                    armChoice.forzed = 1;
                    desc = 'forzed alternation';
                else
                    armChoice.forzed = 0;
                    desc = 'spontaneous alternation';
                end
            end
    
            h = figure;
            hold on
            plot(armChoice.timestamps, armChoice.visitedArm,'color',[.7 .7 .7]);
            scatter(armChoice.timestamps(isnan(armChoice.choice)),...
                armChoice.visitedArm(isnan(armChoice.choice)),100,[.8 .8 .8],'filled');
            scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
                armChoice.visitedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled');
            scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
                armChoice.visitedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled');
            if size(armChoice.delay.timestamps,1)== 2
                armChoice.delay.timestamps = armChoice.delay.timestamps';
            end
            for ii = 1:size(armChoice.delay.timestamps,1)
                fill([armChoice.delay.timestamps(ii,:)'; flip(armChoice.delay.timestamps(ii,:))'],[1 1 1.2 1.2]',...
                [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
            end
            xlabel('seconds'); ylim([-.2 1.2]);
            text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
                desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
                num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
                's'));
            set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});
    
            mkdir('Behavior');
            saveas(h,'Behavior\armChoice.png');
    
            C = strsplit(pwd,'\');
            save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
        else
            warning('DigitalIn format does not match. Was T maze performed? ');
        end  
    else
        % Using manual ttls (when ttls didn't work properly

        tracking = getSessionTracking;
        basepath = basenameFromBasepath();

        if save_original_ttls
            digitalIn = getDigitalIn;
            digitalIn_original = digitalIn;
            save([basepath,'.digitalIn_original.events.mat'],'digitalIn_original');
            digitalIn = [];
            % Copy tracking ttls
            digitalIn.timestampsOn{dlcTtl_channel} = digitalIn_original.timestampsOn{dlcTtl_channel};
            digitalIn.timestampsOff{dlcTtl_channel} = digitalIn_original.timestampsOff{dlcTtl_channel};
            digitalIn.ints{dlcTtl_channel} = digitalIn_original.ints{dlcTtl_channel};
            digitalIn.dur{dlcTtl_channel} = digitalIn_original.dur{dlcTtl_channel};
            digitalIn.intsPeriods{dlcTtl_channel} = digitalIn_original.dur{dlcTtl_channel};
            % Copy fiber ttls
            try
                digitalIn.timestampsOn{fiberTtl_channel} = digitalIn_original.timestampsOn{fiberTtl_channel};
                digitalIn.timestampsOff{fiberTtl_channel} = digitalIn_original.timestampsOff{fiberTtl_channel};
                digitalIn.ints{fiberTtl_channel} = digitalIn_original.ints{fiberTtl_channel};
                digitalIn.dur{fiberTtl_channel} = digitalIn_original.dur{fiberTtl_channel};
                digitalIn.intsperiods{fiberTtl_channel} = digitalIn_original.intsPeriods{fiberTtl_channel};
            catch
            end
        end

        % ROI for left sensor
        if ~isempty(dir('ROI_left.mat'))
            file = dir('ROI_left.mat'); load(file.name);
        else
            disp('Draw ROI for left sensor...');
            h1 = figure;
            imagesc(tracking.avFrame.r);
            axis ij;
            colormap gray;
            roi_left = drawpolygon;
            ROI_left = [roi_left.Position; roi_left.Position(1,:)];
            save('ROI_left.mat','ROI_left');
        end
        
        % ROI for right sensor
        if ~isempty(dir(['ROI_right.mat']))
            file = dir('ROI_right.mat'); load(file.name);
        else
            disp('Draw ROI for right sensor...');
            roi_right = drawpolygon;
            ROI_right = [roi_right.Position; roi_right.Position(1,:)];
            save('ROI_right.mat','ROI_right');
        end

        % ROI for Home Delay
        if ~isempty(dir('ROI_home.mat'))
            file = dir('ROI_home.mat'); load(file.name);
        else
            disp('Draw ROI for Home Delay sensor...');
            roi_home = drawpolygon;
            ROI_home = [roi_home.Position; roi_home.Position(1,:)];
            save('ROI_home.mat','ROI_home');
            close (h1);
        end

        
        h1 = figure;
        hold on;
        imagesc(tracking.avFrame.xSize,tracking.avFrame.ySize,tracking.avFrame.r);
        colormap gray;
        p = plot(ROI_left(:,1)*tracking.convFact, ROI_left(:,2)*tracking.convFact,'r','LineWidth',2);
        p2 = plot(ROI_right(:,1)*tracking.convFact, ROI_right(:,2)*tracking.convFact,'r','LineWidth',2);
        p3 = plot(ROI_home(:,1)*tracking.convFact, ROI_home(:,2)*tracking.convFact,'r','LineWidth',2);
        axis tight;
        axis ij;
        legend(p,'Maze Areas ROI','Location','best');
        xlabel('Normalize/ cm');
        saveas(h1,'Behavior\AreasROI.png');

        % Compute alternations when animal enters the different ROIs
        [in_left,on_left] = inpolygon(tracking.position.x,tracking.position.y,ROI_left(:,1)*tracking.convFact,ROI_left(:,2)*tracking.convFact);
        left_ts = tracking.timestamps(find(diff(in_left) == 1));
        left_ts(2,1:length(left_ts)) = ones(1,length(left_ts))*1;
        

        [in_right,on_right] = inpolygon(tracking.position.x,tracking.position.y,ROI_right(:,1)*tracking.convFact,ROI_right(:,2)*tracking.convFact);
        right_ts = tracking.timestamps(find(diff(in_right) == 1));
        right_ts(2,1:length(right_ts)) = ones(1,length(right_ts))*2;

        [in_home,on_home] = inpolygon(tracking.position.x,tracking.position.y,ROI_home(:,1)*tracking.convFact,ROI_home(:,2)*tracking.convFact);
        home_ts = tracking.timestamps(find(diff(in_home) == 1));
        home_ts(2,1:length(home_ts)) = ones(1,length(home_ts))*3;

        all_ts = [left_ts right_ts home_ts];
        [~,idx] = sort(all_ts(1,:));

        all_ts_sorted = all_ts(:,idx);

        while(all_ts_sorted(2,1)) == 3
            all_ts_sorted(:,1) = [];
        end
        
        ts_area(1) = all_ts_sorted(2,1);
        ts_diff = diff(all_ts_sorted(2,:));

        ts_area = [ts_area ts_diff];
        ts_area_good = find(ts_area ~= 0);

        alternations = all_ts_sorted(:,ts_area_good);

        % Left arm TTls
        digitalIn.timestampsOn{leftArmTtl_channel} = alternations(1,alternations(2,:) == 1)';
        digitalIn.timestampsOff{leftArmTtl_channel} = digitalIn.timestampsOn{leftArmTtl_channel} + 0.3;
        digitalIn.ints{leftArmTtl_channel} = [digitalIn.timestampsOn{leftArmTtl_channel} digitalIn.timestampsOff{leftArmTtl_channel}];
        

        % Right arm TTLs
        digitalIn.timestampsOn{rightArmTtl_channel} = alternations(1,alternations(2,:) == 2)';
        digitalIn.timestampsOff{rightArmTtl_channel} = digitalIn.timestampsOn{rightArmTtl_channel} + 0.3;
        digitalIn.ints{rightArmTtl_channel} = [digitalIn.timestampsOn{rightArmTtl_channel} digitalIn.timestampsOff{rightArmTtl_channel}];

        % Home delay Ttls
        digitalIn.timestampsOn{homeDelayTtl_channel} = alternations(1,alternations(2,:) == 3)';
        digitalIn.timestampsOff{homeDelayTtl_channel} = digitalIn.timestampsOn{homeDelayTtl_channel} + 0.1 ; 
        digitalIn.ints{homeDelayTtl_channel} = [digitalIn.timestampsOn{homeDelayTtl_channel} digitalIn.timestampsOff{homeDelayTtl_channel}];

        

        % score for alternation task
        if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
            armChoice.timestamps = [digitalIn.timestampsOn{leftArmTtl_channel}; digitalIn.timestampsOn{rightArmTtl_channel}]; 
            % 0 is left, 1 is right
            armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{leftArmTtl_channel})); ones(size(digitalIn.timestampsOn{rightArmTtl_channel}))];
            armChoice.delay.timestamps = digitalIn.ints{homeDelayTtl_channel};
            if size(armChoice.visitedArm,1) < size(digitalIn.timestampsOn{homeDelayTtl_channel},1) - 10
                warning('There was problem with one of the sensors! Triying to fix it!')
                leftArmSensor = digitalIn.timestampsOn{leftArmTtl_channel};
                rightArmSensor = digitalIn.timestampsOn{rightArmTtl_channel};
                delaySensor = digitalIn.timestampsOn{homeDelayTtl_channel};
                if size(leftArmSensor,1) < size(rightArmSensor,1) -10
                    warning('Estimating left-arm sensor times...');
                    for ii = 2:length(delaySensor)-1
                        if any(leftArmSensor > delaySensor(ii) & leftArmSensor < delaySensor(ii+1)) | ...
                                any(rightArmSensor > delaySensor(ii) & rightArmSensor < delaySensor(ii+1))
                        else
                            leftArmSensor = [leftArmSensor; mean(delaySensor(ii:ii+1))];
                        end
                    end
                elseif size(rightArmSensor,1) < size(leftArmSensor,1) -10
                    warning('Estimating left-arm sensor times...');
                    for ii = 2:length(delaySensor)-1
                        if any(leftArmSensor > delaySensor(ii) & leftArmSensor < delaySensor(ii+1)) | ...
                                any(rightArmSensor > delaySensor(ii) & rightArmSensor < delaySensor(ii+1))
                        else
                            rightArmSensor = [rightArmSensor; mean(delaySensor(ii:ii+1))];
                        end
                    end
                else error('Not able to guess what is the missing sensor!!');
                end
                armChoice.timestamps = [leftArmSensor; rightArmSensor];
                armChoice.visitedArm = [zeros(size(leftArmSensor)); ones(size(rightArmSensor))];
            end
            
            [armChoice.timestamps, idx] = sort(armChoice.timestamps);
            armChoice.visitedArm = armChoice.visitedArm(idx);
            
            armChoice.delay.dur = nanmean(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1));
            armChoice.delay.durations = round(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1),0);
            armChoice.choice = [NaN; abs(diff(armChoice.visitedArm))]; % 1 is right, 0 is wrong
            armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice) - 1);
            performance = [];
            
            if ~isnan(armChoice.delay.durations(end))
                armChoice.delay.durations(end) = NaN;
            end
            durations = unique(armChoice.delay.durations);
            for ii = 1:length(durations)- 1
                performance = [performance; sum(armChoice.choice(find(armChoice.delay.durations == durations(ii)) + 1)) / length(find(armChoice.delay.durations == durations(ii)))];
            end
            armChoice.delay.performance = performance;
            armChoice.delay.uniqueDurations = durations;
            armChoice.task = task;
            armChoice.expectedArm = [NaN; ~xor(armChoice.visitedArm(2:end), armChoice.choice(2:end))];
    
            if size(digitalIn.timestampsOn,2) >=6
                if digitalIn.timestampsOn{6}>digitalIn.timestampsOff{6}
                    armChoice.forzed = 1;
                    desc = 'forzed alternation';
                else
                    armChoice.forzed = 0;
                    desc = 'spontaneous alternation';
                end
            else
                if armChoice.performance == 1 
                    armChoice.forzed = 1;
                    desc = 'forzed alternation';
                else
                    armChoice.forzed = 0;
                    desc = 'spontaneous alternation';
                end
            end
    
            h = figure;
            hold on
            plot(armChoice.timestamps, armChoice.visitedArm,'color',[.7 .7 .7]);
            scatter(armChoice.timestamps(isnan(armChoice.choice)),...
                armChoice.visitedArm(isnan(armChoice.choice)),100,[.8 .8 .8],'filled');
            scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
                armChoice.visitedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled');
            scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
                armChoice.visitedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled');
            if size(armChoice.delay.timestamps,1)== 2
                armChoice.delay.timestamps = armChoice.delay.timestamps';
            end
            for ii = 1:size(armChoice.delay.timestamps,1)
                fill([armChoice.delay.timestamps(ii,:)'; flip(armChoice.delay.timestamps(ii,:))'],[1 1 1.2 1.2]',...
                [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
            end
            xlabel('seconds'); ylim([-.2 1.2]);
            text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
                desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
                num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
                's'));
            set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});
    
            mkdir('Behavior');
            saveas(h,'Behavior\armChoice.png');
    
            C = strsplit(pwd,'\');
            save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
        else
            warning('DigitalIn format does not match. Was T maze performed? ');
        end  


        




        









        
    end
end

end

%     h1 = figure;
%     subplot(3,1,[1 2])
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     scatter(mazeVirtual(:,1), mazeVirtual(:,2),5,linearizeMazeVirtual);
%     ylabel('cm'); 
%     subplot(3,1,[3])
%     hold on
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMaze)),...
%         vlinMaze,5,vlinMaze);
%     scatter(linspace(0,max(vlinMaze)+50,length(vlinMazeCont)),...
%         vlinMazeCont,5,vlinMazeCont);
%     xlabel('cm');  ylabel('cm');
%     colormap jet
%     saveas(h1,'Behavior\virtualMaze.png');

%     centroidR = [30, 22]; 
%     centroidL = [30, 55];
%     x = behavior.positions.x';
%     y = behavior.positions.y';
%     t = behavior.timestamps;
%     [rangle,rradius] = cart2pol(x-centroidR(1), y - centroidR(2));
%     [langle,lradius] = cart2pol(x-centroidL(1), y - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + .15);
%     langle = wrapTo2Pi(-langle + 5.3658 + .15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial; % Distance traveled in cm, left trial after righ trial
%     
%     h1 = figure;
%     hold on
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     xlabel('cm'); ylabel('cm'); 
%     plot(centroidR(1),centroidR(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     plot(centroidL(1),centroidL(2),'o','MarkerFaceColor',[.8 .4 .3],'MarkerEdgeColor','none');
%     saveas(h1,'Behavior\centroids.png');
%     
%     % 
%     h2 = figure;
%     subplot(3,1,[1 2])
%     scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     prev = 0;
%     linearized = []; arm = [];
%     for ii = 1:size(armChoice.delay.timestamps,1)
%         winTrial = [prev armChoice.delay.timestamps(ii)];
%         prev = armChoice.delay.timestamps(ii);
%         xspam = behavior.timestamps >= winTrial(1) & behavior.timestamps <= winTrial(2);
%         subplot(3,1,[1 2])
%         hold on
%         p = plot(x(xspam),y(xspam),'lineWidth',2);
%         ylabel('cm'); 
%         
%         if armChoice.arm(ii) == 1
%             linearized = [linearized rangle(xspam)];
%             arm = [arm ones(size(rangle(xspam)))];
%         elseif armChoice.arm(ii) == 0
%             linearized = [linearized langle(xspam)];
%             arm = [arm zeros(size(rangle(xspam)))]; 
%         end
%         subplot(3,1,3)
%         hold on
%         plot(t(xspam),linearized(xspam),'lineWidth',2);
%         xlabel('samples');
%         drawnow;
%         pause;    
%         delete(p);
%     end
%     % close(h2);
%       
% %     for ii = 1:size(behavior.timestamps)
% %     end    
% while 1
%     [xg,yg]=ginput(1);
%     fprintf('x: %3.2f, y: %3.2f \n',xg,yg);
%     
%     [rangle,rradius] = cart2pol(xg-centroidR(1), yg - centroidR(2));
%     [langle,lradius] = cart2pol(xg-centroidL(1), yg - centroidL(2));
%     rangle = wrapTo2Pi(rangle - 0.9458 + 0.15);
%     langle = wrapTo2Pi(-langle + 5.3658 + 0.15);
%     
%     rangle = (rad2deg(rangle))/(360/lenghTrial); % +180 to go from 0 to 360
%     langle = (rad2deg(langle))/(360/lenghTrial) + lenghTrial ;
%     fprintf('linearized R: %3.2f \n',rangle);
%     fprintf('linearized L: %3.2f \n\n',langle);
% end
