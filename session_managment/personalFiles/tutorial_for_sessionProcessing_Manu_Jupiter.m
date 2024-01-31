%% Tutorial for session processing_Jupiter

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',1,'analog_optogenetic_channels',[],'promt_hippo_layers',true,'profileType','hippocampus');
close all

% 5% Index session
indexNewSession('copyFiles', true);

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'uLEDiscotheque',...
        'analysis_project_path', adapt_filesep([dropbox_path '\ProjectsOnLine\monoSynBition\data']),'loadLast',false); 
% other projects include InterneuronsLibrary, desVIPnhibition, etc

% PS. uLED sessions
pulses = getAnalogPulses('manualThr',true,'overwrite',true); % 1-index
getDigitalIn;
uLEDPulses = getuLEDPulses('Current',2.7);
processSession('digital_optogenetic_channels',[11 12 13 14 15 16],'analog_optogenetic_channels',[3 4 5 6 7 8],'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false);
% something werid with the pulses... to check!!!
[uLEDResponses] = getuLEDResponse;
NotebookMonoSynBition;
