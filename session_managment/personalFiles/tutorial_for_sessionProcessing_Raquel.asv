%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fPv5', 'Y:\fPv5'},'E:\data\fPv5');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','E:\data\fPv5','analysisPath',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[],'analog_optogenetic_channels',1,'promt_hippo_layers',true,'manual_analog_pulses_threshold',true);


% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'VIPcortex',...
        'analysis_project_path', 'C:\Users\valeg\Dropbox\ProjectsOnLine\desVIPnhibition','loadLast',false);
    
    
    