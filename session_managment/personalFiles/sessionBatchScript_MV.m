%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'U:\fCr1', 'Y:\fCr1'},'X:\data\fCr1');

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

%2% Preprocessing
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

batch_preprocessSession('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

%3% Computing summary
batch_sessionSummary('basepath','E:\data\fCck1','analogChannelsList',[],'digitalChannelsList',1);

batch_sessionSummary('basepath','G:\data\fNkx11','analogChannelsList',65,'digitalChannelsList',[]);


% After Spike sorting, got
processSession('analog_optogenetic_channels',65,'promt_hippo_layers',true,'manual_analog_pulses_threshold',true,'bazler_ttl_channel',1);
 
%% others
editDatFile(pwd,[26460 inf],'option','zeroes');   % remove zeroes
