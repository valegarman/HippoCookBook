%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis computer (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder_temp({'F:\camk12'},'D:\camk12');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
preprocessSession('basepath','D:\camk12\camk12_250508_sess6','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'digitalChannelsList',[1],'sessionSummary',false,'getPos',false);

computeSessionSummary('digitalChannelsList',[1,2]);

uLEDPulses = getuLEDsPulses_legacy();
getuLEDResponse('restrict_to_baseline',false,'uLEDPulses',uLEDPulses);

optogeneticResponses = getOptogeneticResponse('force',true,'digitalChannelsList',[1,2,3,4,5,6,7,8,9,10,11,12,13]);













batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
% ANDREA!!!!!!! add ('LED_threshold',.8) in the function for fSst3!!!!!!!!!!!!
processSession('digital_optogenetic_channels',[1],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

% 5% Revise output using ProcessSession_notebook
edit ProcessSession_notebook.m

% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\Andrea Gallardo\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',true);
    
