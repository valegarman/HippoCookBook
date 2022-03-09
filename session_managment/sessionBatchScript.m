%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

%2% Preprocessing
batch_preprocessSession('basepath','E:\data\fCck1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);


%3% Computing summary
batch_sessionSummary('basepath','E:\data\fCck1','analogChannelsList',1,'digitalChannelsList',1);


preprocessSession('basepath','E:\data\fCck1\fCck1_220215_sess12','cleanArtifacts',({[],1}));

preprocessSession('basepath','D:\Dropbox\DATA\sharedRecordings\NewXmlAnimal\190222\rec1_220219_sess1');

computeSessionSummary('basepath','E:\data\fCck1\fCck1_220217_sess14','analogChannelsList',0,'digitalChannelsList',1);

% Index!
indexNewSession('SWChannel',18,'analogChannelsList',65);
 
%% others
editDatFile(pwd,[0.001 1],'option','zeroes');   % remove zeroes
