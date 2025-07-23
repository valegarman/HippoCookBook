%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis computer (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder_temp({'E:\app2'},'C:\app2');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
% ANEL!!!!!!! add ('LED_threshold',.8) in the function for fSst3!!!!!!!!!!!!
processSession('basepath','C:\wt3\wt3_241120_sess9','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);


% 5% Check out output using ProcessSession_notebook
edit Notebook_for_checking_processSession_results.m

% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\Andrea Gallardo\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',true);
    


% % Ripple dectection
    ExcludeIntervals = [];
    rippleChannel =42; %which is the best ripple channel
    SWChannel = 6;
    eventSpikeThreshold_shanks = [1 2 3 4 5]; % which shanks will be accounted for the spike threshold 
    rippleMasterDetector_threshold = [.50 1]; % [1.5 3.5]
    eventSpikeThreshold = 0; % .5