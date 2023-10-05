function [maze] = getTrials_PMaze(varargin)
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

cd(basepath);
if isempty(tracking)
    tracking = getSessionTracking;
end
% 
% if isempty(armChoice)
%     armChoice = getArmChoice;
% end

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


% Get center of the maze
% center_y = round((max(tracking.position.y) - min(tracking.position.y)) / 2,2);
% center_x = round((max(tracking.position.x) - min(tracking.position.x)) / 2,2);
center_x = tracking.apparatus.centre.x;
center_y = tracking.apparatus.centre.y;
% Getting (0,0) at the center of the maze
tracking.position.y = tracking.position.y - center_y;
tracking.position.x = tracking.position.x - center_x;

% Defining Maze
maze = [];
maze.type = 'PMaze';
maze.radius_in = 29;
maze.pos_y_limits = [-20,30];
maze.pos_x_limits = [-2,2];

maze.polar_rho_limits = [30,75];
maze.polar_theta_limits = [15,120];

maze.boundary{1} = [0,25];
maze.boundary{2} = [0,-20];
maze.boundary{3} = [-10,20];
maze.boundary{4} = [10,20];
maze.boundary{5} = [30,58];
maze.boundary{6} = [20,38];

% Determining polar coordinates
[polar_theta,polar_rho] = cart2pol(tracking.position.x,tracking.position.y);

figure;
plot(tracking.position.x,tracking.position.y,'color',[.5 .5 .5]);
hold on;
figure;
polarplot(polar_theta,polar_rho)

% Changing from polar angle to position along circle by multiplying with the radius (inner radius)
% Negative values are along the left arm, positive values along the right arm
polar_theta = polar_theta*maze.radius_in;


disp('Defining trials for the behavior')
boundary = maze.boundary;

% Determining spatial cross points 
pos1test = find(diff(polar_rho > boundary{1}(2))==-1 & tracking.position.x(1:end-1) < boundary{1}(1)+10 & tracking.position.x(1:end-1) > boundary{1}(1)-10 & tracking.position.y(1:end-1) > 0); % Onset of central arm
pos2test = find(diff(tracking.position.y < boundary{2}(2))==1 & tracking.position.x(1:end-1) < boundary{2}(1)+10 & tracking.position.x(1:end-1) > boundary{2}(1)-10); % End of central arm
pos3test = find(diff(tracking.position.x < boundary{3}(1))==1 & tracking.position.x(1:end-1) < boundary{3}(2)+10 & tracking.position.y(1:end-1) < boundary{3}(2)-10); % Start of left arm
pos4test = find(diff(tracking.position.x > boundary{4}(1))==1 & tracking.position.x(1:end-1) < boundary{4}(2)+10 & tracking.position.y(1:end-1) < boundary{4}(2)-10); % Start of right arm
pos5test = find(diff(polar_theta < boundary{5}(2))== 1 & tracking.position.x(1:end-1) < 5 & polar_rho(1:end-1) > boundary{5}(1)); % End of left rim
pos6test = find(diff(polar_theta < boundary{6}(2))==-1 & tracking.position.x(1:end-1) > 5 & polar_rho(1:end-1) > boundary{6}(1)-5); % End of right rim

pos7home = sort([pos5test;pos6test]);
pos2test(find(diff(pos2test)<tracking.samplingRate)+1) = [];

disp('Plotting position')
figure
plot(tracking.position.x,tracking.position.y,'color',[0.5 0.5 0.5]), hold on
plot([-10,10]+boundary{1}(1),-boundary{1}(2)*[1,1],'r')
plot([-10,10]+boundary{2}(1),boundary{2}(2)*[1,1],'r')
plot([1,1]*boundary{3}(1),-boundary{3}(2)+[-10,10],'r')
plot([1,1]*boundary{4}(1),-boundary{4}(2)+[-10,10],'r')


plot(tracking.position.x(pos1test),tracking.position.y(pos1test),'or') % Onset of central arm
plot(tracking.position.x(pos2test),tracking.position.y(pos2test),'om') % End of central arm
plot(tracking.position.x(pos3test),tracking.position.y(pos3test),'ok') % Start of Left arm
plot(tracking.position.x(pos4test),tracking.position.y(pos4test),'og') % Start of right arm
plot(tracking.position.x(pos5test),tracking.position.y(pos5test),'xc') % Left rim
plot(tracking.position.x(pos6test),tracking.position.y(pos6test),'xy') % Right rim
title('3D position of the animal'), xlabel('X'), ylabel('Y'), zlabel('Z'),axis tight,%view(2)


trials = [];
trials.states.error = []; % Boolean, if a trial is considered an error
trials.states.left_right = []; % Numeric: 1 or 2
trials.stateNames.left_right = {'Left','Right'};
trials.start = 0; % Start time of
trials.end = [];

% Preparing trials matrix
tracking.trials = nan(1,length(tracking.timestamps));
trials_states = zeros(1,length(tracking.timestamps));

i = 0;

% Behavioral scoring
for j = 1:length(pos2test)
    test1 = find(pos1test < pos2test(j));
    test2 = find(pos7home > pos2test(j));
    if ~isempty(test2) & ~isempty(test1)
        if trials.start(end)-pos1test(test1(end)) ~= 0
            i = i+1;
            trials.start(i) = pos1test(test1(end));
            trials.end(i) = pos7home(test2(1));
            tracking.trials(trials.start(i):trials.end(i)) = i;
            trials_states(trials.start(i):trials.end(i)) = 1;
            if sum(ismember(pos5test,trials.end(i)))
                % Left trial
                trials.states.left_right(i) = 1;
                if sum(ismember(pos4test,trials.start(i):trials.end(i)))
                    trials.states.error(i) = true;
                    else
                    trials.states.error(i) = false;
                end
            elseif sum(ismember(pos6test,trials.end(i)))
                % Right trial
                trials.states.left_right(i) = 2;
                if sum(ismember(pos3test,trials.start(i):trials.end(i)))
                    trials.states.error(i) = true;
                else
                    trials.states.error(i) = false;
                end
            else
                trials.states.left_right(i) = 0;
            end
        end
    end
end

% Changning from units of samples to units of time
trials.start = tracking.timestamps(trials.start);
trials.end = tracking.timestamps(trials.end);
trials.nTrials = numel(trials.start);

figure,
subplot(1,2,1)
plot(tracking.timestamps-tracking.timestamps(1),tracking.trials,'.k','linewidth',2), xlabel('Time (sec)'), ylabel('Trials')
subplot(3,2,2)
stairs(tracking.timestamps-tracking.timestamps(1),trials_states,'.-k','linewidth',1), xlabel('Time (sec)'), ylabel('Trial')
subplot(3,2,4)
stairs(trials.states.left_right,'.-b','linewidth',1), xlabel('Trials'), ylabel('Left/Right'), 
yticks(1:numel(trials.stateNames.left_right)), yticklabels(trials.stateNames.left_right)
subplot(3,2,6)
stairs(trials.states.error,'.-k','linewidth',1), xlabel('Trials'), ylabel('Errors')


% SAVEOUTPUT

trials.trials_states = trials_states;
maze.trials = trials;
maze.position.polar_theta = polar_theta;
maze.position.polar_rho = polar_rho;

[~,fbasename,~] = fileparts(pwd);
save([basepath filesep fbasename '.PMaze.Behavior.mat'],'maze');



end

