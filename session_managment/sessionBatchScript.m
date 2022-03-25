%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'U:\fCr1', 'Y:\fCr1'},'X:\data\fCr1');

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

%2% Preprocessing
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

batch_preprocessSession('basepath','G:\data\fPv3','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

%3% Computing summary
batch_sessionSummary('basepath','E:\data\fCck1','analogChannelsList',[],'digitalChannelsList',1);

batch_sessionSummary('basepath','G:\data\fNkx11','analogChannelsList',65,'digitalChannelsList',[]);





preprocessSession('basepath','G:\data\fNkx11\fNkx11_201102_sess13','cleanArtifacts',({65,1}),'analogChannelsList',65);

preprocessSession('basepath','D:\Dropbox\DATA\sharedRecordings\NewXmlAnimal\190222\rec1_220219_sess1');

computeSessionSummary('basepath','G:\data\fNkx11\fNkx11_201113_sess22','analogChannelsList',65,'digitalChannelsList',[]);

% Index!
indexNewSession('analogChannelsList',65,'promt_hippo_layers',true);
 
%% others
editDatFile(pwd,[0.001 1],'option','zeroes');   % remove zeroes
