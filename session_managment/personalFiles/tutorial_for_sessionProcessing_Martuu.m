%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis computer (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'J:\data\fCck1');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',6:16);

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

%% FARE SOLO QUESTO FINO AL 5% 

% 3% CLEAN SESSIONS MANUALLY BY PHY - the one i did using phy2 
  

% 4.1 create session.dat

% values paramters for fCamk10
digitalChannelsList = [6:16];
analogChannelsList =[];
bazler_ttl_channel =5;
anymaze_ttl_channel =[];

session = sessionTemplate(pwd,'showGUI',false);  
session.channels = 1:session.extracellular.nChannels;
session.analysisTags.digital_optogenetic_channels = digitalChannelsList;
session.analysisTags.analog_optogenetic_channels = analogChannelsList;
if ~isempty(bazler_ttl_channel)
    session.analysisTags.bazler_ttl_channel = bazler_ttl_channel;
 end
 if ~isempty(anymaze_ttl_channel)
    session.analysisTags.anymaze_ttl_channel = anymaze_ttl_channel;
 end
 save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
 

% 4.2 Processs individual sessions by by 'processSession'. Example:

processSession('digital_optogenetic_channels',1,'analog_optogenetic_channels',[],'promt_hippo_layers',true);      % digital_optogenetic_channels (0-16) ask manu 
                                                                                                                   % analog_optogenetic_channels - check the folder pulses, if no analog pulses is inside, put []
                                                                                                                   % promt_hippo_layers always put 'true', if analisys is done in cortex put 'false'
% 5% Revise output using ProcessSession_notebook
edit Notebook_for_checking_processSession_results_Martu.m

% 5% Index session
indexNewSession;

%%

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\Andrea Gallardo\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',true);


