function [armChoice] = getCircularMazeArmChoice(varargin)
% Compute T maze performance and gets timestamps
%
% USAGE
%
%   [armChoice, behavior] = getCircularMazeArmChoice(varargin)
%
% INPUTS
% 
% digitalIn                     digitalIn structure with T maze convention:
%                                 1. AnyMaze TTL1,            2. AnyMaze TTL2, 
%                                 3. Left Alternation,        4.Righ Alternation
%                                 5. Home Delay,              6. Is alternation forzed?
%                                 7. 0 is cue left, 1 is cue right (for
%                                       cueSide task)
%
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
%  Pablo Abad 2022. Based on getArmChoice by Manu Valero

%% Defaults and Parms
p = inputParser;
addParameter(p,'digitalIn',[],@ischar);
addParameter(p,'task',[]);
addParameter(p,'force',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'leftArmTtl_channel',3,@isnumeric);
addParameter(p,'rightArmTtl_channel',4,@isnumeric);
addParameter(p,'homeDelayTtl_channel',5,@isnumeric);
addParameter(p,'delayPlot',true,@islogical);

parse(p,varargin{:});
digitalIn = p.Results.digitalIn;
task = p.Results.task;
force = p.Results.force;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;
delayPlot = p.Results.delayPlot;

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

if strcmpi(task,'forced')
    if isempty(digitalIn.timestamps{leftArmTtl_channel})
        try
            file = dir('*_TMaze.txt');
            file_txt = load(file.name);
        catch
            warning('No information about left and right rewards. Quitting...');
        end
    end
    
        
    


elseif strcmpi(task,'alternation')
    % score for alternation task
    if isfield(digitalIn,'timestampsOn') && size(digitalIn.timestampsOn,2)>= 5
        armChoice.timestamps = [digitalIn.timestampsOn{leftArmTtl_channel}'; digitalIn.timestampsOn{rightArmTtl_channel}']; 
        % 0 is left, 1 is right
        armChoice.visitedArm = [zeros(size(digitalIn.timestampsOn{leftArmTtl_channel}))'; ones(size(digitalIn.timestampsOn{rightArmTtl_channel}))'];
        armChoice.delay.timestamps = digitalIn.ints{homeDelayTtl_channel};
%         armChoice.delay.timestamps = digitalIn.timestampsOn{homeDelayTtl_channel};
        armChoice.delay.delay = diff(armChoice.delay.timestamps')';
        
        if size(armChoice.visitedArm,1) < size(digitalIn.timestampsOn{homeDelayTtl_channel},2) - 10
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
        
        try
            armChoice.delay.dur = nanmean(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1));
            armChoice.delay.durations = round(armChoice.delay.timestamps(:,2) - armChoice.delay.timestamps(:,1),0);
            armChoice.delay.durations = [NaN; armChoice.delay.durations];
            armChoice.delay.timestamps = [NaN NaN ;armChoice.delay.timestamps];
%             armChoice.delay.timestamps = [armChoice.delay.timestamps];
            armChoice.choice = [NaN; abs(diff(armChoice.visitedArm))]; % 1 is right, 0 is wrong
            armChoice.performance = nansum(armChoice.choice)/(length(armChoice.choice) - 1);
        
            if length(armChoice.delay.durations) > length(armChoice.choice)
                armChoice.delay.durations(end) = [];
                armChoice.delay.timestamps(end,:) = [];
            end
            durations = unique(armChoice.delay.durations);
            performance = [];
            for ii = 1:length(durations)- 1
                performance = [performance; sum(armChoice.choice(find(armChoice.delay.durations == durations(ii)))) / length(find(armChoice.delay.durations == durations(ii)))];
            end
            armChoice.delay.performance = performance;
            armChoice.delay.uniqueDurations = durations;
        
        catch
            warning('No delay computed...');
        end
                
%         if ~isnan(armChoice.delay.durations(end))
%             armChoice.delay.durations(end) = NaN;
%         end
        
        armChoice.task = task;
        armChoice.expectedArm = [NaN; ~xor(armChoice.visitedArm(2:end), armChoice.choice(2:end))];
%         armChoice.delay.timestamps = armChoice.delay.timestamps(:,1);
        
        armChoice.forzed = 0;
        desc = 'spontanepous alternation';
        
        % General figure ( no including delay)
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
        try
            text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
                desc, ', delay: ',{' '},num2str(round(armChoice.delay.dur,2)), ', # trials: ',{' '},...
                num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
                's'));
        catch
            text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
                desc, ', # trials: ',{' '},...
                num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
                's'));
        end
        set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});

        mkdir('Behavior');
        saveas(h,'Behavior\armChoice.png');

        % Figures including delay
        try
            if length(armChoice.delay.uniqueDurations) > 2
                for jj = 1:length(armChoice.delay.uniqueDurations)-1
                    h = figure;
                    hold on
                    plot(armChoice.timestamps, armChoice.visitedArm,'color',[.7 .7 .7]);
                        scatter(armChoice.timestamps(isnan(armChoice.choice)),...
                    armChoice.visitedArm(isnan(armChoice.choice)),100,[.8 .8 .8],'filled');
                    scatter(armChoice.timestamps(find(armChoice.choice == 1 & armChoice.delay.durations == armChoice.delay.uniqueDurations(jj))),...
                        armChoice.visitedArm(find(armChoice.choice == 1 & armChoice.delay.durations == armChoice.delay.uniqueDurations(jj))),100,[.6 .9 .7],'filled');
                    scatter(armChoice.timestamps(find(armChoice.choice == 0 & armChoice.delay.durations == armChoice.delay.uniqueDurations(jj))),...
                        armChoice.visitedArm(find(armChoice.choice == 0 &  armChoice.delay.durations == armChoice.delay.uniqueDurations(jj))),100,[.9 .6 .7],'filled');
                    ts = find(armChoice.delay.durations == armChoice.delay.uniqueDurations(jj));
                    for ii = 1:length(ts)
                        fill([armChoice.delay.timestamps(ts(ii),:)'; flip(armChoice.delay.timestamps(ts(ii),:))'],[1 1 1.2 1.2]',...
                            [.7 .6 .9],'EdgeColor',[.7 .6 .9],'FaceAlpha',.5)
                    end
                    xlabel('seconds'); ylim([-.2 1.2]);
                    text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.delay.performance(jj),2)),',',{' '},...
                        desc, ', delay: ',{' '},num2str(round(armChoice.delay.uniqueDurations(jj),2)), ', # trials: ',{' '},...
                        num2str(length(armChoice.choice(armChoice.delay.durations == armChoice.delay.uniqueDurations(jj)))),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
                        's'));
                    set(gca,'YTick', [0 1],'YTickLabel',{'Left','Right'});

                    mkdir('Behavior');
                    saveas(h,['Behavior\armChoice_delay_',num2str(armChoice.delay.uniqueDurations(jj)),'.png']);
                end
            end
        end
        
        close all;
        % Save Output
        C = strsplit(pwd,'\');
        save([C{end} '.ArmChoice.Events.mat'], 'armChoice');
        
    else
        warning('DigitalIn format does not match. Was T maze performed? ');
    end  
end

end
