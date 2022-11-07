function [behavior] = getBehaviourYMaze(varargin)
% Creates behaviour for YMaze recordings
%
% USAGE
%
%   [behaviour] = getBehaviourOpenField(varargin)
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
%   Pablo Abad 2022

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
addParameter(p,'useTTLs',false,@islogical);
addParameter(p,'linearize',false,@islogical);

parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;
useTTLs = p.Results.useTTLs;
linearize = p.Results.linearize;


if ~isempty(dir('*YMaze.Behavior.mat')) && ~forceReload 
    disp('YMaze already computed! Loading file.');
    file = dir('*YMaze.Behavior.mat');
    load(file.name);
    return
end

cd(basepath);
if isempty(tracking)
    tracking = anyMazeTracking([],[]);
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

if isfield(tracking,'headposition')
    x_head = tracking.headposition.x;
    y_head = tracking.headposition.y;
end

if isfield(tracking,'tailposition')
    x_tail = tracking.tailposition.x;
    y_tail = tracking.tailposition.y;
end

if isfield(tracking,'zone')
    zone = tracking.zone;
end

% generate events
maps = [];

%     maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
%     maps{ii}(:,2) = linCont(arm==armList(ii));
maps{1}(:,1) = tracking.timestamps;
maps{1}(:,2) = tracking.position.x;
maps{1}(:,3) = tracking.position.y;

%% Get events (animal entering a zone)
try
    % Get entries in a different way
    for ii = 1:length(tracking.zone.bndgbox)
        [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y,tracking.zone.bndgbox{ii}.bndgbox.Vertices(:,1), tracking.zone.bndgbox{ii}.bndgbox.Vertices(:,2));
        a{ii} = diff(in{ii});
        entry.(tracking.zone.bndgbox{ii}.name).sample = find(a{ii} == 1) + 1;
        entry.(tracking.zone.bndgbox{ii}.name).ts = t(entry.(tracking.zone.bndgbox{ii}.name).sample);
        exit.(tracking.zone.bndgbox{ii}.name).sample = find(a{ii} == -1) + 1;
        exit.(tracking.zone.bndgbox{ii}.name).ts = t(exit.(tracking.zone.bndgbox{ii}.name).sample);
    end
       
    h1 = figure;
    plot(tracking.position.x, tracking.position.y, 'Color', [0.5 0.5 0.5])
    hold on;
    axis ij;
    flds = fields(entry);
    plot(tracking.apparatus.bndgbox,'FaceAlpha',0);
    hold on;
    for ii = 1:length(fields(entry)) 
        plot(tracking.zone.bndgbox{ii}.bndgbox,'FaceAlpha',0);
        scatter(tracking.position.x(in{ii}),tracking.position.y(in{ii}),3,'k');
        hold on;
        scatter(tracking.position.x(entry.(flds{ii}).sample),tracking.position.y(entry.(flds{ii}).sample),15,'g');
        scatter(tracking.position.x(exit.(flds{ii}).sample),tracking.position.y(exit.(flds{ii}).sample),15,'r');  
    end    
    xlim(tracking.avFrame.xSize);
    ylim(tracking.avFrame.ySize);
    mkdir('Behavior');
    saveas(h1,'Behavior\Behavior_inZone.png');
    close(h1);
catch
    disp('There are not zones declared. Skipping...');
end

%% Output
behavior.masks = [];
behavior.events = [];
behavior.trials = [];

behavior.timestamps = tracking.timestamps;

behavior.position.lin = nan(length(behavior.timestamps),1);
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
end

behavior.masks.arm = NaN;
behavior.masks.direction = nan(length(behavior.timestamps),1);
behavior.masks.trials = nan(length(behavior.timestamps),1);
behavior.masks.trialsDirection = NaN;


behavior.maps = maps;

try
    behavior.description = tracking.apparatus.name;
catch
    behavior.description = 'YMaze';
end

behavior.events.startPoint = NaN;
behavior.events.rReward = NaN;
behavior.events.lReward = NaN;
behavior.events.startDelay = NaN;
behavior.events.endDelay = NaN;
behavior.events.intersection = NaN;

behavior.trials.startPoint = [NaN NaN];
behavior.trials.endDelay = NaN;
behavior.trials.visitedArm = NaN;
behavior.trials.choice = NaN;
behavior.trials.expectedArm = NaN;

if exist('entry','var')
    try
        behavior.events.entry = entry;
        behavior.events.exit = exit;
    catch
        behavior.events.entry.ts = NaN;
        behavior.events.exit.ts = NaN;
    end
else
    behavior.events.entry.ts = NaN;
    behavior.events.exit.ts = NaN;
end
if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.YMaze.Behavior.mat'], 'behavior');
end

end
