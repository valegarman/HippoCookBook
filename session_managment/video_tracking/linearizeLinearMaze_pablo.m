function [behavior] = linearizeLinearMaze_pablo(varargin)
% Linearize behavior and add relevant events, acording to the euclidean distance to a virtual maze
%
% USAGE
%
%   [behavior] = linearizeLinearMaze(varargin)
%
% INPUTS
% (OPTIONAL)
% basePath            -(default: pwd) basePath for the recording file, in
%                        buzcode format. 
% tracking            - Tracking structure, with a timestamps field and a position field that
%                        contains x (1xC) and y (1xC) subfields. By default, runs LED2Tracking 
%                        to get it.
% digitalIn           - DigitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left arm,          4.Righ arm
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
% behavior.trials.startPoint         Trial epochs, defined as epochs
%                                       between one side to the other side
% behavior.trials.endDelay           Trial epochs, defnied as delays door openings.
% behavior.trials.arm                (1x#trials). Trial's arm (ej 0 left, 1 right)
% 
%   Manu Valero 2020
% Modified by Pablo Abad to depend on session metadata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'editLOI',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'verbose',false,@islogical);
parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
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
%
cd(basepath);
if isempty(tracking)
    tracking = LED2Tracking;
end

if isempty(digitalIn)
    digitalIn = getDigitalIn;
end
if isempty(tracking) || isempty(digitalIn)
    warning('Missing components. No behaviour performed?');
    return
end

% get components
average_frame  = tracking.avFrame.r;        % get average frames
xMaze = tracking.avFrame.xSize;
yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

if isfield(tracking,'zone')
    zone = tracking.zone;
end

cd(basepath); cd ..; upBasepath = pwd; cd(basepath);
if exist([basepath filesep 'virtualMaze.mat'],'file')
    load([basepath filesep 'virtualMaze.mat'],'maze');
elseif exist([upBasepath filesep 'virtualMaze.mat'],'file')
    load([upBasepath filesep 'virtualMaze.mat'],'maze');
    disp('Virtual trajectory from master folder... copying locally...');
    save([basepath filesep 'virtualMaze.mat'],'maze');
else 
    disp('Draw LOI for tracking linearization (one vertex per corner):');
    h0 = figure;
    hold on
    imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 .7]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    caxis([t(1) t(end)]);
    xlim([xMaze]); ylim([yMaze]);
%     xlim([yMaze]); ylim([xMaze]);
    title('Draw a line (hold cursor) following animal trajectory (top to bottom)...','FontWeight','normal');
    maze = drawline;
    maze = [maze.Position];
    editLOI = 'true';
    close(h0);
    save([basepath filesep 'virtualMaze.mat'],'maze');
end

if editLOI
    disp('Edit LOI for tracking linearization:');
    h0 = figure;
    hold on
    imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 .7]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    caxis([t(1) t(end)]);
    xlim([xMaze]); ylim([yMaze]);
%     xlim([yMaze]); ylim([xMaze]);
    title('Move vertex to match trajectory and press Enter...','FontWeight','normal');
    roi = images.roi.Polyline(gca,'Position',maze);
    pause;
    maze = [roi.Position];
    close(h0);
    save([basepath filesep 'virtualMaze.mat'],'maze');
end

% linMazeCont = [0   110]; % left to right
linMazeCont = xMaze; % left to right

% gets steps along the maze
dMaze = diff(maze,1);
dist_vertex = hypot(dMaze(:,1),dMaze(:,2));
cum_dist = [0; cumsum(dist_vertex,1)];
num_points = 5000;
dist_steps = linspace(0, cum_dist(end), num_points);
mazeVirtual = interp1(cum_dist, maze, dist_steps);
vlinMazeCont = interp1(cum_dist, linMazeCont, dist_steps);

disp('Linearizing trajectory...');
for ii = 1:length(x)
    euc = sqrt((mazeVirtual(:,1)-x(ii)).^2 ...
        + (mazeVirtual(:,2)-y(ii)).^2); % euclidean distance between point and virtual maze trajectory
    [~,idEuc] = min(euc);
    linCont(ii) = vlinMazeCont(idEuc);
end

% get trials
try
    cd ..
    upBasepath = pwd;
    session = loadSession();
    cd(basepath)
end
digitalIn = getDigitalIn('all','fs',session.extracellular.sr);
if length(digitalIn.timestampsOn{3}) <  length(digitalIn.timestampsOn{4})- 5 || length(digitalIn.timestampsOn{3}) >  length(digitalIn.timestampsOn{4})+ 5
    warning('Malfunctioning sensor!! Trying to fix');
    if length(digitalIn.timestampsOn{3}) <  length(digitalIn.timestampsOn{4})- 5
        warning('Left sensor is not working properly');
        disp('Find most distant place between consecutive right sensor pulses...');
        rightSensorTimes = digitalIn.timestampsOn{4};
        timestamps = tracking.timestamps;
        rightSensorPositions = interp1(timestamps, linCont', rightSensorTimes);
        for ii = 1:length(rightSensorTimes) - 1
            span = timestamps(timestamps >= rightSensorTimes(ii) & timestamps <= rightSensorTimes(ii+1));
            [~, idx] = max(abs(mean(rightSensorPositions(ii:ii+1)) - linCont(timestamps >= rightSensorTimes(ii) & timestamps <= rightSensorTimes(ii+1))));            
            leftSensorTimes(ii,1) =  span(idx);
        end
        leftSensorPositions = interp1(timestamps, linCont', leftSensorTimes);
    else
        warning('Right sensor is not working properly');
        disp('Find most distant place between consecutive right sensor pulses...');
        leftSensorTimes = digitalIn.timestampsOn{3};
        timestamps = tracking.timestamps;
        leftSensorPositions = interp1(timestamps, linCont', leftSensorTimes);
        for ii = 1:length(leftSensorTimes) - 1
            span = timestamps(timestamps >= leftSensorTimes(ii) & timestamps <= leftSensorTimes(ii+1));
            [~, idx] = max(abs(mean(leftSensorPositions(ii:ii+1)) - linCont(timestamps >= leftSensorTimes(ii) & timestamps <= leftSensorTimes(ii+1))));            
            rightSensorTimes(ii,1) =  span(idx);
        end
    end
else
    rightSensorTimes = digitalIn.timestampsOn{4};
    leftSensorTimes = digitalIn.timestampsOn{3};
end

rightSensorTimes = rightSensorTimes';
leftSensorTimes = leftSensorTimes';

armChoice.timestamps = [rightSensorTimes; leftSensorTimes];
armChoice.visitedArm = [ones(size(rightSensorTimes)); zeros(size(leftSensorTimes))];
[armChoice.timestamps, idx] = sort(armChoice.timestamps);
armChoice.visitedArm = armChoice.visitedArm(idx);

% armChoice.timestamps = armChoice.timestamps';
% armChoice.visitedArm = armChoice.visitedArm';

h2 = figure;
subplot(3,1,[1 2])
hold on
% scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
plot(mazeVirtual(:,1), mazeVirtual(:,2),'k-');
colormap parula
colorTraj = jet(size(armChoice.timestamps,1));
prev = 0; 
direction = ones(size(linCont)); % left to right is 1
winTrial = []; 
arm = ones(size(linCont));

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
        plot(t(xspam),linCont(xspam),'color',colorTraj(ii,:),'lineWidth',1);
        % plot(t(xspam),lin(xspam),'g','lineWidth',2);
        xlabel('seconds'); ylabel('cm');

        if armChoice.visitedArm(ii) == 0
            arm(xspam) = 0;
        end
        plot(tracking.timestamps(xspam([end end])),[0 max(y)],'k');

        subplot(3,1,[1 2])
        if verbose    
            drawnow;
            pause(0.01); 
            delete(p);
        end
    end
end
try
    for ii = 1:length(zone.bndgbox)
        plot(zone.bndgbox{ii}.bndgbox,'FaceAlpha',0);
    end
catch
end
plot(mazeVirtual(:,1), mazeVirtual(:,2),'k-');

% interpolate events
rReward = rightSensorTimes;
for ii = 1:length(rReward)
    [~,idx] = min(abs(rReward(ii) - t));
    p3 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .1],'MarkerEdgeColor','k');
end
lReward = leftSensorTimes;
for ii = 1:length(lReward)
    [~,idx] = min(abs(lReward(ii) - t));
    p4 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.1 .5 .8],'MarkerEdgeColor','k');
end

legend([p3 p4],'rReward', 'lReward');

saveas(h2,'Behavior\linearizeTrajectory.png');

if ~verbose
    close(h2);
end

% generate events
maps = [];
armList = unique(arm);
for ii = 1:length(armList)
    maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
    maps{ii}(:,2) = linCont(arm==armList(ii));
end

trials0 = [armChoice.timestamps [armChoice.timestamps(2:end); t(end)]]; % trials defined as epochs between sensors crossing

% generate trials mask
trialMask = nan(size(tracking.timestamps));
for ii = 1:length(trials0)
    posTrials = find(tracking.timestamps >= trials0(ii,1) & tracking.timestamps <= trials0(ii,2));
    trialMask(posTrials) = ii;
    trialMaskDirection(ii) = round(mean(arm(posTrials)));
end

% populate behavior
behavior.timestamps = tracking.timestamps;

behavior.position.lin = linCont';
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
end

behavior.masks.arm = arm';
behavior.masks.direction = arm';
behavior.masks.trials = trialMask;
behavior.masks.trialsDirection = trialMaskDirection;

behavior.maps = maps;

try
    behavior.description = tracking.apparatus.name;
catch
    behavior.description = 'linearMaze';
end

behavior.events.startPoint = NaN;
behavior.events.rReward = rReward;
behavior.events.lReward = lReward;
behavior.events.startDelay = NaN;
behavior.events.endDelay = NaN;
behavior.events.intersection = NaN;
behavior.events.entry.ts = NaN;
behavior.events.exit.ts = NaN;

behavior.trials.startPoint = trials0;
behavior.trials.endDelay = NaN;
behavior.trials.visitedArm = armChoice.visitedArm;
behavior.trials.choice = NaN;
behavior.trials.expectedArm = NaN;


if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Linearized.Behavior.mat'], 'behavior');
end


end
