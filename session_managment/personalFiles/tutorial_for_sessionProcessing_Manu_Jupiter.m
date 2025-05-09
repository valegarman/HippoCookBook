%% Tutorial for session processing_Jupiter

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','Y:\unindexedSubjects\fCamk10','analysisPath','C:\Data\fCamk10','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',6:16,'bazler_ttl_channel',5,'skipStimulationPeriods',false);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% <OPTIONAL> 
preprocessSession('basepath',pwd,'analogChannelsList',[],'spikeSort',true,'getPos',false, 'cleanArtifacts',false,...
                    'medianSubstr',true,'sessionSummary',true,'digitalChannelsList',[1:16]);


% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[1],'analog_optogenetic_channels',[],'promt_hippo_layers',true,'profileType','hippocampus');
close all

% 5% Index session
indexNewSession('copyFiles', false);

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'Bibliocampus',...
        'analysis_project_path', adapt_filesep([onedrive_path 'NeuralComputationLab\ActiveProjects\interneuronsLibrary\data']),'loadLast',false); 
% other projects include Bibliocampus, desVIPnhibition, etc

% PS. uLED sessions
pulses = getAnalogPulses('manualThr',true); % 1-index
getDigitalIn;
uLEDPulses = getuLEDPulses('Current',2.7);
processSession('digital_optogenetic_channels',[1:16],'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false);
% something werid with the pulses... to check!!!
[uLEDResponses] = getuLEDResponse;
NotebookMonoSynBition;


% 
getOptogeneticResponse('digitalChannelsList', 9, 'force', true)
