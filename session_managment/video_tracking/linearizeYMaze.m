function [behavior] = linearizeYMaze(varargin)
% Linearize behavior and add relevant events, acording to the euclidean distance to a virtual maze
%
% USAGE
%
%   [behavior] = linearizeYMaze(varargin)
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
%   Manu Valero 2019

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
    tracking = LED2Tracking;
end

if isempty(armChoice)
    armChoice = getYMazeArmChoice;
end

if isempty(digitalIn)
    digitalIn = getDigitalIn;
end

if isempty(tracking) || isempty(digitalIn)
    warning('Missing components. No behaviour performed?');
    return
end

if isfield(tracking,'zone')
    zone = tracking.zone;
end

% get components
average_frame  = tracking.avFrame.r;        % get average frames
xMaze = tracking.avFrame.xSize;
yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

cd(basepath); cd ..; upBasepath = pwd; cd ..; up_upBasepath = pwd; cd(basepath);
if exist([basepath filesep 'virtualMaze.mat'],'file')
    load([basepath filesep 'virtualMaze.mat'],'maze');
elseif exist([upBasepath filesep 'virtualMaze.mat'],'file')
    load([upBasepath filesep 'virtualMaze.mat'],'maze');
    disp('Virtual trajectory from session folder... copying locally...');
    save([basepath filesep 'virtualMaze.mat'],'maze');
elseif exist([up_upBasepath filesep 'virtualMaze.mat'],'file')
    load([up_upBasepath filesep 'virtualMaze.mat'],'maze');
    disp('Virtual trajectory from master folder... copying locally...');
    save([basepath filesep 'virtualMaze.mat'],'maze');
else 
    disp('Draw LOI for tracking linearization (Center first to right):');
    h0 = figure;
    hold on
    axis ij;
    imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 .7]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    caxis([t(1) t(end)]);
    xlim([xMaze]); ylim([yMaze]);
    title('Draw a polyline following animal trajectory (Starting from center to right)...','FontWeight','normal');
    maze = drawpolyline;
    maze = [maze.Position; maze.Position(1,:)];
    editLOI = 'true';
    close(h0);
    save([basepath filesep 'virtualMaze.mat'],'maze');
end

if editLOI
    disp('Edit LOI for tracking linearization:');
    h0 = figure;
    hold on
    axis ij;
    imagesc(xMaze, yMaze,average_frame); colormap gray; caxis([0 .7]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    caxis([t(1) t(end)]);
    xlim([xMaze]); ylim([yMaze]);
    title('Move vertex to match trajectory and press Enter...','FontWeight','normal');
    roi = images.roi.Polyline(gca,'Position',maze);
    pause;
    maze = [roi.Position];
    close(h0);
    save([basepath filesep 'virtualMaze.mat'],'maze');
end

linMazeCont =    [0   45  ... % right turn
                 0       ... % right arm to center
                 45      ... % left turn
                 0      ... % left arm to center
                 45      ... % stem arm
                 0      ... % stem arm to center
                 ];         
             
% gets steps along the maze
dMaze = diff(maze,1);
dist_vertex = hypot(dMaze(:,1),dMaze(:,2));
cum_dist = [0; cumsum(dist_vertex,1)];
num_points = 5000;
dist_steps = linspace(0, cum_dist(end), num_points);
mazeVirtual = interp1(cum_dist, maze, dist_steps);
vlinMazeCont = interp1(cum_dist, linMazeCont, dist_steps);

% correct center linearization

vlinMazeCont(vlinMazeCont>=linMazeCont(5)) = ...
    vlinMazeCont(vlinMazeCont>=linMazeCont(5)) - linMazeCont(5); 
    
disp('Linearizing trajectory...');
for ii = 1:length(x)
    euc = sqrt((mazeVirtual(:,1)-x(ii)).^2 ...
        + (mazeVirtual(:,2)-y(ii)).^2); % euclidean distance between point and virtual maze trajectory
    [~,idEuc] = min(euc);
    linCont(ii) = vlinMazeCont(idEuc);
end

% add path distane to correct for stam lineaization when go left
% stem = find(linCont<=linMazeCont(2));   % timestamp stem
% rarm = find(linCont>linMazeCont(3) & linCont < linMazeCont(4));        % timestamp r arm
% larm = find(linCont>linMazeCont(7) & linCont< linMazeCont(8));        % timestamp l arm
% for ii = 1:length(stem)
%     rt =(rarm - stem(ii)); rt(rt<=0) = NaN;
%     lt =(larm - stem(ii)); lt(lt<=0) = NaN;
%     if min(lt) < min(rt) % if in stem and go left
%         linCont(stem(ii)) = linCont(stem(ii)) + 170;
%     end
% end

% check trajectories, generate intersection point timestamp, check left
% and right reward...

%% Figure 
h2 = figure;
subplot(3,1,[1 2])
hold on
axis ij;
% scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
plot(mazeVirtual(:,1), mazeVirtual(:,2),'k-');
colormap parula
colorTraj = jet(size(t,1));
prev = 0; 
arm = zeros(size(linCont));
traj = zeros(size(linCont));
winTrial = [];

for ii = 1:size(armChoice.timestamps,1)-1
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
        xlabel('seconds');

        if armChoice.visitedArm(ii) == 2
            arm(xspam) = 2;
        elseif armChoice.visitedArm(ii) == 1
            arm(xspam) = 1;
        end
        
        [~,idxInt] = min(abs(linCont(xspam)-0)); % find point closer to 45 in lin
        intersection(ii) = t(xspam(idxInt));
        sampleIntersection(ii) = xspam(idxInt);

        plot([intersection(ii) intersection(ii)],[0 linMazeCont(2)],'k');

        subplot(3,1,[1 2])
        p1 = plot(x(sampleIntersection(ii)),y(sampleIntersection(ii)),'o'...
            ,'MarkerFaceColor','w','MarkerEdgeColor','k');

        if verbose    
            drawnow;
            pause(0.01); 
            delete(p);
        end
    end
end
plot(mazeVirtual(:,1), mazeVirtual(:,2),'k-');
    
% interpolate events
stemArm = armChoice.timestamps(armChoice.visitedArm==0);
for ii = 1:length(stemArm)
    [~,idx] = min(abs(stemArm(ii) - t));
    p2 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .1],'MarkerEdgeColor','k');
end
rightArm = armChoice.timestamps(armChoice.visitedArm==1);
for ii = 1:length(rightArm)
    [~,idx] = min(abs(rightArm(ii) - t));
    p3 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.1 .5 .8],'MarkerEdgeColor','k');
end
leftArm = armChoice.timestamps(armChoice.visitedArm==2);
for ii = 1:length(leftArm)
    [~,idx] = min(abs(leftArm(ii) - t));
    p4 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .8],'MarkerEdgeColor','k');
end
legend([p1 p2 p3 p4],'Intersection', 'StemArm', 'rightArm', 'leftArm');
    
saveas(h2,'Behavior\linearizeTrajectoryYMaze.png');

if ~verbose
    close(h2);
end
    

%% Figure

h2 = figure;
subplot(3,1,[1 2])
hold on
axis ij;
% scatter(x,y,3,[.8 .8 .8],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
plot(mazeVirtual(:,1), mazeVirtual(:,2),'k-');
colormap parula
colorTraj = jet(size(t,1));
prev = 0; 
arm = zeros(size(linCont));
traj = zeros(size(linCont));
winTrial = [];

colormap jet
colorTraj0 = jet(size(t,1));
colorTraj1 = jet(size(t,1))*2 / 255;
colorTraj2 = jet(size(t,1))*3 / 255;
colorTraj3 = jet(size(t,1))*4 / 255;
colorTraj4 = jet(size(t,1))*5 / 255;
colorTraj5 = jet(size(t,1))*6 / 255;


for ii = 1:size(armChoice.timestamps,1)-1
    winTrial = [prev armChoice.timestamps(ii)];
    prev = armChoice.timestamps(ii);
    xspam = find(tracking.timestamps >= winTrial(1) & tracking.timestamps <= winTrial(2));
    
    if ~isempty(xspam)
       subplot(3,1,[1 2])
       hold on;
       if armChoice.trajectory(ii) == 0 % stem to right
            traj(xspam) = 0;
            scatter(x(xspam),y(xspam),10,colorTraj0(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
       elseif armChoice.trajectory(ii) == 1 % right to stem
            traj(xspam) = 1;
            scatter(x(xspam),y(xspam),10,colorTraj1(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
        elseif armChoice.trajectory(ii) == 2 % right to left
            traj(xspam) = 2;
            scatter(x(xspam),y(xspam),10,colorTraj2(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
        elseif armChoice.trajectory(ii) == 3 % left to right
            traj(xspam) = 3;
            scatter(x(xspam),y(xspam),10,colorTraj3(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
       elseif armChoice.trajectory(ii) == 4 % stem to left
           traj(xspam) = 4;
           scatter(x(xspam),y(xspam),10,colorTraj4(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
       elseif armChoice.trajectory(ii) == 5 % left to right
           traj(xspam) = 5;
           scatter(x(xspam),y(xspam),10,colorTraj5(ii,:),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
       elseif isnan(armChoice.trajectory(ii))
           traj(xspam) = NaN;
           scatter(x(xspam),y(xspam),10,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
       end
       
       ylabel('cm');
       subplot(3,1,3)
       hold on;
       
        
        
    end



end

%% generate events
maps = [];
armList = unique(arm);
for ii = 1:length(armList)
    maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
    maps{ii}(:,2) = linCont(arm==armList(ii));
end

maps_whole{1}(:,1) = tracking.timestamps;
maps_whole{1}(:,2) = linCont;

trials0 = [intersection' [(intersection(2:end))'; t(end)]]; % trials defined as epochs between 0 position crossings

% generate trials mask
trialMask = nan(size(tracking.timestamps));
for ii = 1:length(trials0)
    posTrials = find(tracking.timestamps >= trials0(ii,1) & tracking.timestamps <= trials0(ii,2));
    trialMask(posTrials) = ii;
end

% populate behavior
behavior.timestamps = tracking.timestamps;

behavior.position.lin = linCont';
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

behavior.masks.arm = arm';
behavior.masks.trials = trialMask;
behavior.masks.direction = arm';
behavior.masks.trialsDirection = NaN;

behavior.maps = maps;
behavior.maps_whole = maps;

behavior.description = armChoice.task;

behavior.events.startPoint = NaN;
behavior.events.rReward = NaN;
behavior.events.lReward = NaN;
behavior.events.startDelay = NaN;
behavior.events.endDelay = NaN;

behavior.events.stemArm = stemArm;
behavior.events.rightArm = rightArm;
behavior.events.leftArm = leftArm;
behavior.events.intersection = intersection';
behavior.events.entry.ts = NaN;
behavior.events.exit.ts = NaN;

behavior.trials.startPoint = trials0;
behavior.trials.endDelay = NaN;

behavior.trials.visitedArm = armChoice.visitedArm;
behavior.trials.choice = armChoice.choice;
behavior.trials.expectedArm = NaN;

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
else
    behavior.zone = [];
end
if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Linearized.Behavior.mat'], 'behavior');
end

end
