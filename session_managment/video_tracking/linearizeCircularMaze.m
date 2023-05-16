function [behavior] = linearizeCircularMaze(varargin)
% Linearize behavior and add relevant events, acording to the euclidean distance to a virtual maze
%
% USAGE
%
%   [behavior] = linearizeArmChoice(varargin)
%
% INPUTS
% (OPTIONAL)
% basePath            -(default: pwd) basePath for the recording file, in
%                        buzcode format. 
% tracking            - Tracking structure, with a timestamps field and a position field that
%                        contains x (1xC) and y (1xC) subfields. By default, runs LED2Tracking 
%                        to get it.
% armChoice           - Event structure after getArmchoice function             
% digitalIn           - DigitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left Alternation,  4.Righ Alternation
%                                 5. Home Delay,        6. Is alternation forzed?
% editLOI             - Edit loaded Line of interest (LOI). 
% saveMat             - Default true
% forceReload         - Default false
% verbose             - Default true
% 
% OUTPUT
%                     - Behavior structure with the following fields updated:
% 
% behavior.timestamps                Total behavioral timestamps
% behavior.position.lin              Linearized position in cm
% behavior.position.x                X coordinates of tracking, in cm/norm
% behavior.position.y                Y coordinates, in cm/norm 
% behavior.masks.arm                 Code for map maze arms (ej, 0 is left, 1 is arm)
% behavior.maps                      Cell array as [time position], one cell/map
% behavior.description               
% behavior.events
% behavior.trials.startPoint         Trial epochs, defined as epochs between 0 position crossings
% behavior.trials.endDelay           Trial epochs, defnied as delays door openings.
% behavior.trials.arm                (1x#trials). Trial's arm (ej 0 left, 1 right)
% behavior.trial.choice              (1x#trials). 0 is wrong, 1 is right.
% behavior.trial.side                (1x#trial ). In CueSide maze, maps the right side
% 
%   Pablo Abad 2022. Based on linearizeArmChoice by Manuel Valero 2019.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'armChoice',[],@isstruct);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'editLOI',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'verbose',false,@islogical);
parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
armChoice = p.Results.armChoice;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;


if ~isempty(dir('*Linearized.Behavior.mat')) && ~forceReload 
    disp('Linearization already computed! Loading file.');
    file = dir('*Linearized.Behavior.mat');
    load(file.name);
    return
end

% maze = getTrials_PMaze();
maze = getTrials_PMaze2('forceReload',true);

% Arm choice

cd(basepath);
if isempty(tracking)
    tracking = LED2Tracking;
    center_x = tracking.apparatus.centre.x;
    center_y = tracking.apparatus.centre.y;
    tracking.position.x = tracking.position.x - center_x;
    tracking.position.y = tracking.position.y - center_y;
end

if isempty(armChoice)
    armChoice = getArmChoice;
end

if isempty(digitalIn)
    digitalIn = getDigitalIn;
end

if isempty(tracking)  || isempty(digitalIn)
    warning('Missing components. No behaviour performed?');
    return
end

% if isempty(tracking) || isempty(armChoice) || isempty(digitalIn)
%     warning('Missing components. No behaviour performed?');
%     return
% end

close all;

% Defining Arm and rim as states


behavior.states.arm_rim = nan(1,length(tracking.timestamps))'; % Numeric: 1, 2 or nan (outside)
behavior.stateNames.arm_rim = {'Arm','Rim'};

% idx_arm = maze.position.polar_rho < maze.polar_rho_limits(1) & tracking.position.x > maze.pos_x_limits(1)-3 & tracking.position.x < maze.pos_x_limits(2);
idx_arm = tracking.position.x > maze.pos_x_limits(1)-3 & tracking.position.x < maze.pos_x_limits(2);
behavior.states.arm_rim(idx_arm) = 1; % arm positions

% Adding homeport to central arm
idx = find(isnan(behavior.states.arm_rim) & tracking.position.y >= maze.pos_y_limits(2));
pos_linearized(idx) = tracking.position.y(idx)+ maze.pos_y_limits(2);
behavior.states.arm_rim(idx) = 1;

idx_rim = ~idx_arm & maze.position.polar_rho > maze.polar_rho_limits(1) & maze.position.polar_rho < maze.polar_rho_limits(2)+10 &...   
    ((maze.position.polar_theta < 0 & maze.position.polar_theta < -maze.boundary{5}(2)) | ...
    (maze.position.polar_theta > 0 & maze.position.polar_theta > maze.boundary{6}(2)));
behavior.states.arm_rim(idx_rim) = 2;

figure;
plot(tracking.position.x,tracking.position.y,'color',[.5 .5 .5]);
hold on;
scatter(tracking.position.x(idx_arm), tracking.position.y(idx_arm),'k');
scatter(tracking.position.x(idx), tracking.position.y(idx),'k');
scatter(tracking.position.x(idx_rim), tracking.position.y(idx_rim),'r');

% First the central arm is linearized
pos_linearized = nan(size(tracking.timestamps));
pos_linearized(behavior.states.arm_rim==1) = tracking.position.y(behavior.states.arm_rim==1)+maze.pos_y_limits(2);
% pos_linearized(behavior.states.arm_rim==1) = tracking.position.y(behavior.states.arm_rim==1);

% Finally the right return side-arm is linearized.
pos_linearized(behavior.states.arm_rim==2 & maze.position.polar_theta > 0) = (maze.position.polar_theta(behavior.states.arm_rim==2 & maze.position.polar_theta > 0)-90) ;
figure;
plot(pos_linearized)

% Next the left(?) return side-arm
boundary = diff(maze.pos_y_limits);
pos_linearized(behavior.states.arm_rim==2 & maze.position.polar_theta < 0) = -(maze.position.polar_theta(behavior.states.arm_rim==2 & maze.position.polar_theta < 0) + 90);

figure;
plot(pos_linearized)

% check trajectories, generate intersection point timestamp, check left
% and right reward...

% get components
average_frame  = tracking.avFrame.r;        % get average frames
xMaze = tracking.avFrame.xSize;
yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

h2 = figure;
subplot(3,1,[1 2])
hold on
colormap parula
colorTraj = jet(size(armChoice.timestamps,1));
prev = 0; arm = ones(size(pos_linearized));
winTrial = [];
for ii = 1:size(armChoice.timestamps,1)
    winTrial = [prev armChoice.timestamps(ii)];
    prev = armChoice.timestamps(ii);
    xspam = find(tracking.timestamps >= winTrial(1) & tracking.timestamps <= winTrial(2));
    if ~isempty(xspam)
        subplot(3,1,[1 2])
        hold on
        scatter(x(xspam),y(xspam),10,colorTraj(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
        if verbose
            p = scatter(x(xspam),y(xspam),10,colorTraj(ii,:),'filled');
        end
        ylabel('cm'); 

        subplot(3,1,3)
        hold on
        plot(t(xspam),pos_linearized(xspam),'color',colorTraj(ii,:),'lineWidth',1);
        % plot(t(xspam),lin(xspam),'g','lineWidth',2);
        xlabel('seconds');

        if armChoice.visitedArm(ii) == 0
            arm(xspam) = 0;
        end
        
        [~,idxInt] = min(abs(pos_linearized(xspam)-0)); % find point closer to 0 in lin  (intersection)
        intersection(ii) = t(xspam(idxInt));
        sampleIntersection(ii) = xspam(idxInt);

        [~,idxInt] = min(abs(pos_linearized(xspam)-60)); % find point closer to 0 in lin (homeDelay)
        homeCage(ii) = t(xspam(idxInt));
        samplehomeCage(ii) = xspam(idxInt);
        
        if sampleIntersection(ii) == samplehomeCage(ii)
            homeCage(ii) = NaN;
            samplehomeCage(ii) = NaN;
        end

        plot([intersection(ii) intersection(ii)],[0 60],'k');
        if ~isnan(homeCage(ii))
            plot([homeCage(ii) homeCage(ii)],[0 60],'r');
        end

        subplot(3,1,[1 2])
        p1 = plot(x(sampleIntersection(ii)),y(sampleIntersection(ii)),'o'...
            ,'MarkerFaceColor','w','MarkerEdgeColor','k');
        if ~isnan(samplehomeCage(ii))
            p2 = plot(x(samplehomeCage(ii)),y(samplehomeCage(ii)),'o'...
                ,'MarkerFaceColor','k','MarkerEdgeColor','w');
        end

        if verbose    
            drawnow;
            pause(0.01); 
            delete(p);
        end
    end
end

% interpolate events
% rReward = digitalIn.timestampsOn{4};
rReward = armChoice.timestamps(armChoice.visitedArm==1);
for ii = 1:length(rReward)
    [~,idx] = min(abs(rReward(ii) - t));
    p3 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .1],'MarkerEdgeColor','k');
end
% lReward = digitalIn.timestampsOn{3};
lReward = armChoice.timestamps(armChoice.visitedArm==0);
for ii = 1:length(lReward)
    [~,idx] = min(abs(lReward(ii) - t));
    p4 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.1 .5 .8],'MarkerEdgeColor','k');
end
endDelay = armChoice.delay.timestamps(:,2);
for ii = 1:length(endDelay)
    [~,idx] = min(abs(endDelay(ii) - t));
    p5 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .8],'MarkerEdgeColor','k');
end
startDelay = armChoice.delay.timestamps(:,1);
for ii = 1:length(startDelay)
    [~,idx] = min(abs(startDelay(ii) - t));
    p6 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.5 .8 .5],'MarkerEdgeColor','k');
end
legend([p1 p2 p3 p4 p5 p6],'Inters', 'HomeCage', 'rReward', 'lReward', 'endDelay', 'startDelay');
    
saveas(h2,'Behavior\linearizeTrajectory.png');

if ~verbose
    close(h2);
end
    
% generate events
maps = [];
armList = unique(arm);
for ii = 1:length(armList)
    maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
    maps{ii}(:,2) = pos_linearized(arm==armList(ii));
end

endDelay(1) = t(1);
endDelay = endDelay';
startDelay(1) = t(1);
startDelay = startDelay';

% endDelay = [t(1) endDelay];
% startDelay = [t(1) startDelay];

trials0 = [homeCage' [(homeCage(2:end))'; t(end)]]; % trials defined as epochs between 0 position crossings
trialsDelay = [endDelay' [(endDelay(2:end))'; t(end)]]; % trials defined as epochs between end delays positions

% generate trials mask
trialMask = nan(size(tracking.timestamps));
for ii = 1:length(trials0)
    posTrials = find(tracking.timestamps >= trials0(ii,1) & tracking.timestamps <= trials0(ii,2));
    trialMask(posTrials) = ii;
end

% populate behavior
behavior.timestamps = tracking.timestamps;

behavior.position.lin = pos_linearized;
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

behavior.masks.arm = arm';
behavior.masks.trials = trialMask;
behavior.masks.direction = arm';
behavior.masks.trialsDirection = NaN;

behavior.maps = maps;

behavior.description = 'PMaze';
behavior.task = armChoice.task;

behavior.events.startPoint = homeCage';
behavior.events.rReward = rReward;
behavior.events.lReward = lReward;
behavior.events.startDelay = startDelay';
behavior.events.endDelay = endDelay';
behavior.events.intersection = intersection';
behavior.events.entry.ts = NaN;
behavior.events.exit.ts = NaN;

behavior.trials.startPoint = trials0;
behavior.trials.endDelay = trialsDelay;

behavior.trials.visitedArm = armChoice.visitedArm;
behavior.trials.choice = armChoice.choice;
behavior.trials.expectedArm = armChoice.expectedArm;

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Linearized.Behavior.mat'], 'behavior');
end

end