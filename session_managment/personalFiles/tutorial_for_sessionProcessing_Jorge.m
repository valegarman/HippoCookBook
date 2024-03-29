%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fVip5', 'U:\fVip5'},'G:\data\fVip5');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','G:\data\fVip5','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% 3% CLEAN SESSIONS MANUALLY BY PHY


% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'bazler_ttl_channel',10,'promt_hippo_layers',true);


% 5% Index session 
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\manu\Documents','loadLast',false);
    
%% Fix digitalIn not recorded in initial recording
% 1. Run getDigitalIn in subfolder with digitalIn
getDigitalIn;
% 2. Copy/paste the generated amplifier.DigitalIn.events to session folder.
% 3. Run
load('amplifier.DigitalIn.events.mat');
load([basenameFromBasepath(pwd) '.MergePoints.events.mat']);

digitalIn.ints{1} = digitalIn.ints{1} + MergePoints.timestamps(2,1);
digitalIn.intsPeriods{1} = digitalIn.intsPeriods{1} + MergePoints.timestamps(2,1);
digitalIn.timestampsOff{1} = digitalIn.timestampsOff{1} + MergePoints.timestamps(2,1);
digitalIn.timestampsOn{1} = digitalIn.timestampsOn{1} + MergePoints.timestamps(2,1);

save([basenameFromBasepath(pwd) '.DigitalIn.events.mat'],'digitalIn');


    