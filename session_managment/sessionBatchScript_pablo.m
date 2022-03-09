%1% Transfer files and organize session's folder
updateExpFolder({'X:\fCr1', 'Y:\fCr1'},'Z:\data\fCr1');
updateExpFolder({'X:\fCr2'},'Z:\data\fCr2');
preprocessSession('basepath','Z:\data\fCck3\fCck3_220228_sess7','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','Z:\data\fCr1\fCr1_220225_sess1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
batch_preprocessSession('basepath','Z:\data\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% Cck3
updateExpFolder({'V:\data\fCck3', 'Y:\fCck3'},'Z:\data\fCck3');
batch_preprocessSession('basepath','Z:\data\fCck3','analysisPath','C:\data\fCck3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% batch_sessionSummary('basepath','C:\data\fCck3','analogChannelsList',0,'digitalChannelsList',1);

% Cr1
updateExpFolder({'X:\fCr1', 'Y:\fCr1'},'Z:\data\fCr1');
batch_preprocessSession('basepath','Z:\data\fCr1','analysisPath','C:\data\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);


batch_sessionSummary('basepath','E:\data\fCck1','analogChannelsList',1,'digitalChannelsList',1);

% Prueba
preprocessSession('basepath','Z:\data\fCck3\fCck3_220222_sess3','analysisPath','C:\data\fCck3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
