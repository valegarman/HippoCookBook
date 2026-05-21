function [tracking] = getBpod(varargin)
%
% Gets Bpod (Sanworks) information (events and states)   
%
% USAGE
%
%   [Bpod] = getSessionBpod(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file
%
% OUTPUT
%   Bpod.behaviour output structure, with fields:
%       - performance.
%
%   HISTORY:
%     - Pablo Abad 2026 (NeuralComputationLab).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',false,@islogical)

parse(p,varargin{:});

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;


%% Deal with inputs
if ~isempty(dir([basepath filesep '*Bpod.Behavior.mat'])) && forceReload
    disp('Bpod already detected! Loading file.');
    file = dir([basepath filesep '*Bpod.Behavior.mat']);
    load(file.name);
    return
end

%% Load MATLAB file
file = dir('*CentralPokeReward*mat');
load(file.name);

Bpod = [];

Bpod.info.eventNames = SessionData.Info.FSMsetup.StateMachineInfo.EventNames;
Bpod.info.date = SessionData.Info.SessionDate;
Bpod.info.startTime_UTC = SessionData.Info.SessionStartTime_UTC;
Bpod.info.startTime_Matlab = SessionData.Info.SessionStartTime_MATLAB;

flds = fields(SessionData.SettingsFile.GUI);
for ii = 1:length(flds)
    Bpod.settings.(flds{ii}) = SessionData.SettingsFile.GUI.(flds{ii});
end
Bpod.settings.state_names = SessionData.RawData.OriginalStateNamesByNumber{1};

Bpod.numTrials = SessionData.nTrials;

for ii = 1:SessionData.nTrials
    Bpod.trials.num_poke(ii) = length(find(SessionData.RawData.OriginalEventData{ii} == 112));
    if Bpod.trials.num_poke(ii) > 0
        Bpod.trials.latency(ii) = SessionData.RawData.OriginalEventTimestamps{ii}(find(SessionData.RawData.OriginalEventData{ii} == 112,1));
    else
        Bpod.trials.latency(ii) = NaN;
    end
    Bpod.trials.trial_duration(ii) = SessionData.TrialEndTimestamp(ii) - SessionData.TrialStartTimestamp(ii);
    Bpod.trials.start(ii) = SessionData.TrialStartTimestamp(ii);
    Bpod.trials.end(ii) = SessionData.TrialEndTimestamp(ii);
end

Bpod.ITI = SessionData.ITI;
Bpod.response = SessionData.Response;
Bpod.reward = SessionData.Reward;

if saveMat
    filename = [];
end










