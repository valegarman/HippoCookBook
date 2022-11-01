%% Tutorial for session processing

% 0% In case data was recorded in spike2 (or other formats), transform to
% dat file
spike2toDat;

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'F:\fCamk10','Y:\fCamk10'},'J:\fCamk10');
updateExpFolder({'V:\data\fId4','Y:\fId4'},'J:\fId4');
updateExpFolder({'V:\data\fVip4'},'J:\fVip4');


% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','J:\fCamk10','cleanArtifacts',({[],[5 6 7  8 9 10  11 12 13  14 15 16]}),'analogChannelsList',[],'digitalChannelsList',1);
batch_preprocessSession('basepath','J:\fVip4','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);


% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','J:\fVip4','analogChannelsList',[],'digitalChannelsList',1);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[],'analog_optogenetic_channels',1,'promt_hippo_layers',true,'manual_analog_pulses_threshold',false);

processSession('digital_optogenetic_channels',[],'analog_optogenetic_channels',1,'promt_hippo_layers',true,'excludeAnalysis',{'3'});



% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\valeg\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',false);
    