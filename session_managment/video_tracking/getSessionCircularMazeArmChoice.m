function [sessionCircularMazeArmChoice] = getSessionCircularMazeArmChoice(varargin)
% Compute and get session arm choice over all session subfolders
%
% USAGE
%
%   [sessionArmChoice] = getCircularMazeArmChoice(varargin)
%
% INPUTS
% basePath                      (default: pwd) basePath for the recording file, 
%                                    in buzcode format:
% task                          'alternation' and 'cudeSide' 
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
% forceRun                      Try to run armChoice even if there is no
%                                   tracking files on the directories (default false).
%
% OUTPUT
%       - armChoice.behaviour.(subSessionFolder) output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.arm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.ints          Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.delay.timestamps    Delay timestamps, in seconds
% armChoice.choice              Performance vector, 1 is right choice, 0 is
%                                   wrong. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
% armChoice.task                'alternation' and 'cudeSide'  
%
%   Pablo Abad 2022. Based on getSessionArmChoice by Manuel Valero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'task',[],@ischar);
addParameter(p,'force',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceRun',false,@islogical);
addParameter(p,'leftArmTtl_channel',3,@isnumeric);
addParameter(p,'rightArmTtl_channel',4,@isnumeric);
addParameter(p,'homeDelayTtl_channel',5,@isnumeric);

parse(p,varargin{:});
task = p.Results.task;
forceReload = p.Results.force;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
forceRun = p.Results.forceRun;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.SessionArmChoice.Events.mat'])) || forceReload
    disp('Session arm choice already detected! Loading file.');
    file = dir([basepath filesep '*.SessionArmChoice.Events.mat']);
    load(file.name);
    return
end

%% Find subfolder recordings
cd(basepath);

session = loadSession(basepath);

count = 1;
for ii = 1:length(session.epochs)
    if ~isempty(dir([basepath filesep session.epochs{ii}.name filesep '*tracking.csv'])) || forceRun 
        cd([basepath filesep session.epochs{ii}.name]);
        fprintf('Computing CircularMaze Arm Choice in %s folder \n',session.epochs{ii}.name);
        tracking = getSessionTracking();
        if strcmpi(tracking.apparatus.name,'TMaze')
            sessionArmChoice.(session.epochs{ii}.name)= getCircularMazeArmChoice('verbose',verbose,'task',task,'leftArmTtl_channel',leftArmTtl_channel,...
                'rightArmTtl_channel',rightArmTtl_channel,'homeDelayTtl_channel',homeDelayTtl_channel);
            trackFolder(count) = ii; 
            count = count + 1;
        end
    end
end
cd(basepath);

efields = fieldnames(sessionArmChoice);
try tracking = getSessionTracking;
    if size(tracking.events.subSessions,1) == size(efields,1)
        disp('Correctiong timestamps for session recording...');
        for ii = 1:size(efields,1)
            preRec = tracking.events.subSessions(ii,1);
            sessionArmChoice.(efields{ii}).timestamps = ...
                sessionArmChoice.(efields{ii}).timestamps + preRec;
            sessionArmChoice.(efields{ii}).delay.timestamps = ...
                sessionArmChoice.(efields{ii}).delay.timestamps + preRec;
        end
    else
        warning('Number of behavioral recordings do not match!')
    end
catch 
    warning('No available tracking!');
end

try [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
catch
    sessionInfo.FileName = split(pwd,filesep); sessionInfo.FileName = sessionInfo.FileName{end};
end
if saveMat
    save([basepath filesep sessionInfo.FileName '.SessionArmChoice.Events.mat'],'sessionArmChoice');
end

end
