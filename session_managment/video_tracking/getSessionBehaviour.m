
function [behaviour] = getSessionBehaviour(varargin)
% Compute beahviour over all session subfolders for OpenField recordings.
%
% USAGE
%   [behavior] = getSessionBehaviour(varargin)

% INPUTS
% basepath                      (default: pwd) basepath for the recording file, 
%                                    in buzcode format.
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
% maze                          OpenField or YMaze. 
%
% OUTPUT
%
%   Pablo Abad 2021. Based on getSessionLinearize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'verbose',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'maze','OpenField',@ischar)
parse(p,varargin{:});
forceReload = p.Results.forceReload;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
maze = p.Results.maze;

%% Get session metada
session = loadSession(basepath);
%% Deal with inputs
if ~isempty(dir([basepath filesep session.general.name '.Behavior.mat'])) || forceReload
    disp('Behavior already detected! Loading file.');
    file =dir([basepath filesep session.general.name '.Behavior.mat']);
    load(file.name);
    return
end

% if isempty(maze)
%     fileTarget = dir('*SessionArmChoice*');
%     if isempty(fileTarget)
%         maze = 'linearMaze';
%     else
%         maze = 'tMaze';
%     end
% end

% I need to automatize this part in order to get the YMaze or the OpenField
% position. (YMaze should have a SessionArmChoice maybe ¿?¿?¿)

%% Find subfolder recordings
cd(basepath);
C = strsplit(session.general.name,'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 1;
for ii = 1:size(sess,1)
    if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing behaviour in %s folder \n',sess(ii).name);
        
        if ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
            file = dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']);
            load(file.name);
        end
        
        if strcmpi(tracking.apparatus.name,'yMaze') || strcmpi(tracking.apparatus.name,'YMaze Apparatus')
            behaviorTemp.(sess(ii).name)= getBehaviorYMaze; %% I need to create this function for the Ymaze paradigm
        elseif strcmpi(tracking.apparatus.name,'Open Field') || strcmpi(tracking.apparatus.name,'OpenField') || strcmpi(tracking.apparatus.name,'Social Interaction')
            behaviorTemp.(sess(ii).name)= getBehaviourOpenField;
        end
        trackFolder(count) = ii; 
        count = count + 1;
    end
end
cd(basepath);

efields = fieldnames(behaviorTemp);
tracking = getSessionTracking;
x = []; y = []; timestamps = []; 
x_head = []; y_head = [];
x_tail = []; y_tail = [];
lin = []; armMask = []; trialMask = []; recMask = [];
startPoint = []; rReward = []; lReward = []; startDelay = []; endDelay = []; intersection = [];
startPointTrials = []; endDelayTrials = []; visitedArm = []; choice = []; expectedArm = [];
recordings = []; recordingsTrial = []; description = [];
direction = []; trialsDirection = []; events = [];
subSessionMask = [];

try
    for ii = 1:size(efields,1)
        if ~isfield(behaviorTemp.(efields{ii}).masks, 'direction')
            behaviorTemp.(efields{ii}).masks.direction = behaviorTemp.(efields{ii}).masks.arm;
            behaviorTemp.(efields{ii}).masks.trialsDirection = behaviorTemp.(efields{ii}).masks.trials;
        end    
    end
catch
end

% SubSessionsMask
for ii=1:size(efields,1)
    behaviorTemp.(efields{ii}).events.subSessionMask = (tracking.events.subSessionsMask == ii)*ii;
    behaviorTemp.(efields{ii}).events.subSessionMask(behaviorTemp.(efields{ii}).events.subSessionMask == 0) = [];
end

if size(tracking.events.subSessions,1) == size(efields,1)
    disp('Correctiong timestamps for session recording...');
    for ii = 1:size(efields,1)
        preRec = tracking.events.subSessions(ii,1);
        timestamps = [timestamps; behaviorTemp.(efields{ii}).timestamps + preRec];
%         timestamps{ii} = behaviorTemp.(efields{ii}).timestamps + preRec;
        subSessionMask = [subSessionMask;behaviorTemp.(efields{ii}).events.subSessionMask];
        x = [x; behaviorTemp.(efields{ii}).position.x];
        y = [y; behaviorTemp.(efields{ii}).position.y];
%         x{ii} = behaviorTemp.(efields{ii}).position.x;
%         y{ii} = behaviorTemp.(efields{ii}.position.y);
        if isfield(behaviorTemp.(efields{ii}),'headposition')
            x_head = [x_head; behaviorTemp.(efields{ii}).headposition.x];
            y_head = [y_head; behaviorTemp.(efields{ii}).headposition.y];
%             x_head{ii} = behaviorTemp.(efields{ii}).headposition.x];
%             y_head{ii} = behaviorTemp.(efields{ii}).headposition.y];
        end
        if isfield(behaviorTemp.(efields{ii}),'tailposition')
            x_tail = [x_tail; behaviorTemp.(efields{ii}).tailposition.x];
            y_tail = [y_tail; behaviorTemp.(efields{ii}).tailposition.y];
        end
        if isfield(behaviorTemp.(efields{ii}),'zone')
            zone{ii} = behaviorTemp.(efields{ii}).zone;
        end
            
        events{ii} = behaviorTemp.(efields{ii}).events;
        
        lin = [lin; behaviorTemp.(efields{ii}).position.lin];
        try
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
        catch
        end
        description{ii} = behaviorTemp.(efields{ii}).description;
    end
else
    warning('Number of behavioral recordings do not match!')
end


% generate maps, one for each arm and for each recording
maps = [];
count = 1;
for ii = 1:length(efields)
%         maps{count}(:,1) = timestamps(direction==directionList(jj) & recMask==ii);
%         maps{count}(:,2) = lin(direction==directionList(jj) & recMask==ii);
    maps{count}(:,1) = timestamps(subSessionMask == ii);
    maps{count}(:,2) = x(subSessionMask == ii);
    maps{count}(:,3) = y(subSessionMask == ii);
    count = count + 1;
end

% populate behavior
behaviour.timestamps = timestamps;

behaviour.position.lin = lin;
behaviour.position.x = x;
behaviour.position.y = y;

if exist('x_head','var') && ~isempty(x_head)
    behaviour.headposition.x = x_head;
    behaviour.headposition.y = y_head;
end

if exist('x_tail','var') && ~isempty(x_tail)
    behaviour.tailposition.x = x_tail;
    behaviour.tailposition.y = y_tail;
end

if exist('zone','var') && ~isempty(zone)
    behaviour.zone = zone;
end

behaviour.maps = maps;

behaviour.description = description;

behaviour.events = events;

try
    behaviour.masks.arm = armMask;
    behaviour.masks.trials = trialMask;
    behaviour.masks.direction = direction;
    behaviour.masks.trialsDirection = trialsDirection;
    behaviour.masks.recording = recMask;

    behaviour.events.startPoint = startPoint;
    behaviour.events.rReward = rReward;
    behaviour.events.lReward = lReward;
    behaviour.events.startDelay = startDelay;
    behaviour.events.endDelay = endDelay;
    behaviour.events.intersection = intersection;
    behaviour.events.recordings = recordings;

    behaviour.trials.startPoint = startPointTrials;
    behaviour.trials.endDelay = endDelayTrials;
    behaviour.trials.visitedArm = visitedArm;
    behaviour.trials.choice = choice;
    behaviour.trials.expectedArm = expectedArm;
    behaviour.trials.recordings = recordingsTrial;
catch
end

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Behavior.mat'], 'behaviour');
end




end

