
function [behavior] = getSessionBehavior(varargin)
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
addParameter(p,'maze',[],@ischar)
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
            behaviorTemp.(sess(ii).name)= getBehaviourYMaze; %% I need to create this function for the Ymaze paradigm
        elseif strcmpi(tracking.apparatus.name,'Open Field') || strcmpi(tracking.apparatus.name,'OpenField') || strcmpi(tracking.apparatus.name,'Social Interaction')
            behaviorTemp.(sess(ii).name) = getBehaviourOpenField;
        elseif strcmpi(tracking.apparatus.name,'Linear Track  N-S')
            behaviorTemp.(sess(ii).name) = linearizeLinearMaze_pablo('verbose',verbose);
        elseif strcmpi(tracking.apparatus.name,'TMaze')
            behaviorTemp.(sess(ii).name) = linearizeArmChoice('verbose',verbose);
        end
        behaviorFolder(count) = ii; 
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
        if isfield(behaviorTemp.(efields{ii}),'zone')
            zone{ii} = behaviorTemp.(efields{ii}).zone;
        end
            
%         events{ii} = behaviorTemp.(efields{ii}).events;
        
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
%         description{ii} = behaviorTemp.(efields{ii}).description;
    end
else
    warning('Number of behavioral recordings do not match!')
end

% generate maps, one for each arm and for each recording
maps = [];
directionList = unique(direction);
directionList(find(isnan(directionList))) = [];
count = 1;
for ii = 1:length(efields)
    if strcmpi(behaviorTemp.(efields{ii}).description,'Linear Track  N-S') && ~any(find(isnan(behaviorTemp.(efields{ii}).position.lin)))
        for jj = 1:length(directionList)
            maps{count}(:,1) = timestamps(direction == directionList(jj) & recMask == ii);
            maps{count}(:,2) = lin(direction == directionList(jj) & recMask==ii);
            description{count} = behaviorTemp.(efields{ii}).description;
            zone{count} = behaviorTemp.(efields{ii}).zone;
            avFrame{count} = tracking.avFrame{ii};
            count = count+1;
        end
    elseif strcmpi(behaviorTemp.(efields{ii}).description,'Open Field') && any(isnan(behaviorTemp.(efields{ii}).position.lin))
        maps{count}(:,1) = timestamps(subSessionMask == ii);
        maps{count}(:,2) = x(subSessionMask == ii);
        maps{count}(:,3) = y(subSessionMask == ii);
        description{count} = behaviorTemp.(efields{ii}).description;
        zone{count} = behaviorTemp.(efields{ii}).zone;
        avFrame{count} = tracking.avFrame{ii};
        count = count + 1;
    else
        disp('Error while running getSessionBehavior. Quitting...');
        return;
    end
end

% populate behavior
behavior.timestamps = timestamps;

behavior.position.lin = lin;
behavior.position.x = x;
behavior.position.y = y;

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
end

behavior.maps = maps;

behavior.description = description;
behavior.avFrame = avFrame;

for ii = 1:length(events)
    if isfield(events{ii},'entry') && isfield(events{ii},'exit')
        flds = fields(events{ii});
        rmv = find(ismember(flds,{'exit','entry'}) == 0);
        flds{rmv} = [];
        preRec = tracking.events.subSessions(ii,1);
        for jj = 1:length(flds)
            if ~isempty(flds{jj})
                events{ii}.(flds{jj}).leftArm.ts = events{ii}.(flds{jj}).leftArm.ts + preRec;
                events{ii}.(flds{jj}).rightArm.ts = events{ii}.(flds{jj}).rightArm.ts + preRec;
                events{ii}.(flds{jj}).stemArm.ts = events{ii}.(flds{jj}).stemArm.ts + preRec;
                events{ii}.(flds{jj}).centerArm.ts = events{ii}.(flds{jj}).centerArm.ts + preRec;
            end
        end
    end
end


% behavior.events = events;

try
    behavior.masks.arm = armMask;
    behavior.masks.trials = trialMask;
    behavior.masks.direction = direction;
    behavior.masks.trialsDirection = trialsDirection;
    behavior.masks.recording = recMask;

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
catch
end

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Behavior.mat'], 'behavior');
end




end

