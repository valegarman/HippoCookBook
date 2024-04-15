
function [behavior] = getSessionLinearize(varargin)
% Compute behavior linearization over all session subfolders
%
% USAGE
%
%   [behavior] = getSessionLinearizeTmaze(varargin)
%
% INPUTS
% basePath                      (default: pwd) basePath for the recording file, 
%                                    in buzcode format:
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
% maze                          tMaze or linearMaze. If SessionArmChoice is
%                                   exists, it runs Tmaze, otherwise
%                                   linearMaze.
%
% OUTPUT
%
%   Manuel Valero 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'verbose',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'maze',[],@ischar)
addParameter(p,'leftTtl',3,@isnumeric)
addParameter(p,'rightTtl',4,@isnumeric)

parse(p,varargin{:});
forceReload = p.Results.forceReload;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
maze = p.Results.maze;
leftTtl = p.Results.leftTtl;
rightTtl = p.Results.rightTtl;


%% Deal with inputs
filename = split(pwd,filesep); filename = filename{end};
if ~isempty(dir([basepath filesep filename '.Behavior.mat'])) && ~forceReload
    disp('Behavior linearization already detected! Loading file.');
    file =dir([basepath filesep filename '.Behavior.mat']);
    load(file.name);
    return
end

if isempty(maze)
    fileTarget = dir('*SessionArmChoice*');
    if isempty(fileTarget)
        maze = 'linearMaze';
    else
        maze = 'tMaze';
    end
end
    
%% Find subfolder recordings
cd(basepath);
sess = basenameFromBasepath(pwd);
% [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
C = strsplit(basenameFromBasepath(pwd),'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 1;
for ii = 1:size(sess,1)
     if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing linearization in %s folder \n',sess(ii).name);
        if strcmpi(maze,'tMaze')
            behaviorTemp.(sess(ii).name)= linearizeArmChoice('verbose',true);
        elseif strcmpi(maze,'linearMaze')
             disp('sono qui');
            behaviorTemp.(sess(ii).name)= linearizeLinearMaze('verbose',verbose, 'leftTtl',leftTtl ,'rightTtl',rightTtl);
      end
        trackFolder(count) = ii; 
        count = count + 1;
     else
         behaviorTemp=[];
     end
end
cd(basepath);

efields = fieldnames(behaviorTemp);
tracking = getSessionTracking;
x = []; y = []; timestamps = []; lin = []; armMask = []; trialMask = []; recMask = [];
startPoint = []; rReward = []; lReward = []; startDelay = []; endDelay = []; intersection = [];
startPointTrials = []; endDelayTrials = []; visitedArm = []; choice = []; expectedArm = [];
recordings = []; recordingsTrial = []; description = [];
direction = []; trialsDirection = [];

for ii = 1:size(efields,1)
    if ~isfield(behaviorTemp.(efields{ii}).masks, 'direction')
        behaviorTemp.(efields{ii}).masks.direction = behaviorTemp.(efields{ii}).masks.arm;
        behaviorTemp.(efields{ii}).masks.trialsDirection = behaviorTemp.(efields{ii}).masks.trials;
    end
end

if size(tracking.events.subSessions,1) == size(efields,1)
    disp('Correctiong timestamps for session recording...');
    for ii = 1:size(efields,1)
        preRec = tracking.events.subSessions(ii,1);
        timestamps = [timestamps; behaviorTemp.(efields{ii}).timestamps + preRec];
        x = [x; behaviorTemp.(efields{ii}).position.x];
        y = [y; behaviorTemp.(efields{ii}).position.y];
        lin = [lin; behaviorTemp.(efields{ii}).position.lin];
        armMask = [armMask; behaviorTemp.(efields{ii}).masks.arm];
        trialMask = [trialMask; behaviorTemp.(efields{ii}).masks.trials];
        recMask = [recMask; ii * ones(size(behaviorTemp.(efields{ii}).masks.trials))];
        startPoint = [startPoint; behaviorTemp.(efields{ii}).events.startPoint + preRec];
        rReward = [rReward; behaviorTemp.(efields{ii}).events.rReward + preRec];
        lReward = [lReward; behaviorTemp.(efields{ii}).events.lReward + preRec];
        startDelay = [startDelay; behaviorTemp.(efields{ii}).events.startDelay + preRec];
        endDelay = [endDelay; behaviorTemp.(efields{ii}).events.endDelay + preRec];
        intersection = [intersection; behaviorTemp.(efields{ii}).events.intersection + preRec];
        recordings = [recordings; preRec];
        startPointTrials = [startPointTrials; behaviorTemp.(efields{ii}).trials.startPoint + preRec];
        endDelayTrials = [endDelayTrials; behaviorTemp.(efields{ii}).trials.endDelay + preRec];
        visitedArm = [visitedArm; behaviorTemp.(efields{ii}).trials.visitedArm];
        choice = [choice; behaviorTemp.(efields{ii}).trials.choice];
        expectedArm = [expectedArm; behaviorTemp.(efields{ii}).trials.expectedArm];
        direction = [direction; behaviorTemp.(efields{ii}).masks.direction];
        trialsDirection = [trialsDirection; behaviorTemp.(efields{ii}).masks.trialsDirection'];
        recordingsTrial = [recordingsTrial; ii*ones(size(behaviorTemp.(efields{ii}).trials.startPoint,1),1)];
        description{ii} = behaviorTemp.(efields{ii}).description;
    end
else
    warning('Number of behavioral recordings do not match!')
end

% generate maps, one for each arm and for each recording
maps = [];
directionList = unique(direction);
if size(direction,1) ~= 1
    direction = direction';
end
if size(recMask,1) ~= 1
    recMask = recMask';
end
count = 1;
for ii = 1:length(efields)
    for jj = 1:length(directionList)
        maps{count}(:,1) = timestamps(direction==directionList(jj) & recMask==ii);
        maps{count}(:,2) = lin(direction==directionList(jj) & recMask==ii);
        count = count + 1;
    end
end

% populate behavior
if size(timestamps,2) ~= 1
    timestamps = timestamps';
end
behavior.timestamps = timestamps;

behavior.position.lin = lin;
behavior.position.x = x;
behavior.position.y = y;

behavior.masks.arm = armMask;
behavior.masks.trials = trialMask;
behavior.masks.direction = direction;
behavior.masks.trialsDirection = trialsDirection;
behavior.masks.recording = recMask;

behavior.maps = maps;

behavior.description = description;

behavior.events.startPoint = startPoint;
behavior.events.rReward = rReward;
behavior.events.lReward = lReward;
behavior.events.startDelay = startDelay;
behavior.events.endDelay = endDelay;
behavior.events.intersection = intersection;
behavior.events.recordings = recordings;

behavior.trials.startPoint = startPointTrials;
behavior.trials.endDelay = endDelayTrials;
behavior.trials.visitedArm = visitedArm;
behavior.trials.choice = choice;
behavior.trials.expectedArm = expectedArm;
behavior.trials.recordings = recordingsTrial;

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Behavior.mat'], 'behavior');
end

end