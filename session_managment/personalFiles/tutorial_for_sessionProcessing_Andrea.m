%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis computer (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder_temp({'G:\camk13'},'D:\camk13');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
% preprocessSession('basepath','D:\camk13\camk13_250613_sess2\','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[1, 2]}),'digitalChannelsList',[1 2],'sessionSummary',false,'getPos',false, 'medianSubstr', [1:3 4 11:12 14:18 21:23 27 30 32]);
preprocessSession('basepath','Y:\unindexedSubjects\cancer2\cancer2_250703_sess2','analysisPath','E:\','exclude_shanks',[],'cleanArtifacts',[],'digitalChannelsList',[1:8],'sessionSummary',true,'getPos',false, 'medianSubstr', true);

computeSessionSummary('digitalChannelsList',[1,2]);

uLEDPulses = getuLEDsPulses_legacy('uLEDs_ttl',[1 2]);
getuLEDResponse('restrict_to_baseline',false,'uLEDPulses',uLEDPulses);

optogeneticResponses = getOptogeneticResponse('force',true,'digitalChannelsList',[1,2,3,4,5,6,7,8,9,10,11,12,13]);



batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
% ANDREA!!!!!!! add ('LED_threshold',.8) in the function for fSst3!!!!!!!!!!!!
processSession('digital_optogenetic_channels',[1,2,3,4],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

% 5% Revise output using ProcessSession_notebook
edit Notebook_for_checking_processSession_results

% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\Andrea Gallardo\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',true);
    
