%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'V:\data\fPv5','Y:\fPv5'},'E:\data\fPv5');

%2% Preprocessing
% batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

batch_preprocessSession('basepath','E:\data\fPv5','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

%3% Computing summary
batch_sessionSummary('basepath','E:\data\fPv5','analogChannelsList',[],'digitalChannelsList',1);



preprocessSession('basepath','G:\data\fNkx11\fNkx11_201102_sess13','cleanArtifacts',({65,1}),'analogChannelsList',65);

preprocessSession('basepath','D:\Dropbox\DATA\sharedRecordings\NewXmlAnimal\190222\rec1_220219_sess1');

computeSessionSummary('basepath',pwd,'analogChannelsList',[],'digitalChannelsList',1);

% Index!
indexNewSession('analogChannelsList',65,'promt_hippo_layers',true,'manual_analog_pulses_threshold',true,'bazler_ttl_channel',1);
 
%% others
editDatFile(pwd,[0.001 1],'option','zeroes');   % remove zeroes
