%% sessionBatchScript

%%1% Transfer files and organize session's folder
updateExpFolder({'Here type source path!!'},'Here type target path!!'); % for example updateExpFolder({'U:\fCr1', 'Y:\fCr1'},'X:\data\fCr1');

%%2% Preprocessing
batch_preprocessSession('basepath','Here write path containing all recorded session from one subject!!','analysisPath','Here analysis path, ideally and SSD!!',...
    'analogChannelsList','Did you record any analog channel??','digitalChannelsList', 'Any digital channel??');
% for example: batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

%%3% Computing summary
batch_sessionSummary('basepath','Here write path containing all recorded session from one subject!!',...
    'analogChannelsList','Did you record any analog channel??','digitalChannelsList','Any digital channel??');
% for example: batch_sessionSummary('basepath','G:\data\fNkx11','analogChannelsList',65,'digitalChannelsList',[]);

%%4% Index!
%run indexNewSession() for each session
