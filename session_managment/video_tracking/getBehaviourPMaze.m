function [behavior] = getBehaviourPMaze(varargin)
% Creates behaviour for PMaze recordings in 2D (not linearize)
%
% USAGE
%
%   [behaviour] = getBehaviourPMaze(varargin)
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

parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;
useTTLs = p.Results.useTTLs;

if ~isempty(dir('*PMaze.Behavior.mat')) && ~forceReload 
    disp('PMaze already computed! Loading file.');
    file = dir('*PMaze.Behavior.mat');
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
% average_frame  = tracking.avFrame.r;        % get average frames
% xMaze = tracking.avFrame.xSize;
% yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;








% generate events
maps = [];
maps{1}(:,1) = tracking.timestamps;
maps{1}(:,2) = tracking.position.x;
maps{1}(:,3) = tracking.position.y;



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
behavior.maps_whole = NaN;

% behaviour.description = 'Open Field';
try
    behavior.description = tracking.apparatus.name;
catch
    behavior.description = 'PMaze';
end

behavior.events.startPoint = NaN;
behavior.events.rReward = NaN;
behavior.events.lReward = NaN;
behavior.events.startDelay = NaN;
behavior.events.endDelay = NaN;
behavior.events.intersection = NaN;
behavior.events.stemArm = NaN;
behavior.events.rightArm = NaN;
behavior.events.leftArm = NaN;

behavior.trials.startPoint = [NaN NaN];
behavior.trials.endDelay = NaN;
behavior.trials.visitedArm = NaN;
behavior.trials.choice = NaN;
behavior.trials.expectedArm = NaN;

if exist('entry','var')
    try
        behavior.events.entry = entry;
        behavior.events.exit = exit;
%         behavior.events.entry_sample = entry_sample;
%         behavior.events.exit_sample = exit_sample;
%         behavior.events.entry_ts = entry_ts;
%         behavior.events.exit_ts = exit_ts;
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
    save([C{end} '.PMaze.Behavior.mat'], 'behavior');
end

end
