%1% Transfer files and organize session's folder
% updateExpFolder({'V:\data\fCck3', 'Y:\fCck3'},'Z:\data\fCr1');

updateExpFolder({'X:\fCr2'},'Z:\data\fCr2');
preprocessSession('basepath','Z:\data\fCck3\fCck3_220228_sess7','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','Z:\data\fCr1\fCr1_220225_sess1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
batch_preprocessSession('basepath','Z:\data\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% Cck3
updateExpFolder({'V:\data\fCck3', 'Y:\fCck3'},'Z:\data\fCck3');
batch_preprocessSession('basepath','Z:\data\fCck3','analysisPath','C:\data\fCck3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
batch_sessionSummary('basepath','Z:\data\fCck3','analogChannelsList',[],'digitalChannelsList',1);

preprocessSession('basepath','C:\data\fCck3\fCck3_220314_sess17','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','C:\data\fCck3\fCck3_220315_sess18','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','C:\data\fCck3\fCck3_220316_sess19','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','C:\data\fCck3\fCck3_220317_sess20','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

% batch_sessionSummary('basepath','C:\data\fCck3','analogChannelsList',0,'digitalChannelsList',1);

% Cr1
updateExpFolder({'X:\fCr1', 'Y:\fCr1'},'Z:\data\fCr1');
batch_preprocessSession('basepath','Z:\data\fCr1','analysisPath','C:\data\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);


batch_sessionSummary('basepath','E:\data\fCck1','analogChannelsList',1,'digitalChannelsList',1);
% Prueba con fNkx11 median subtsraction a ver si es mas facil el clustering
preprocessSession('basepath','D:\data\fNkx11\fNkx11_201104_sess15','analogChannelsList',65,'digitalChannelsList',[],'manualThr',true);
preprocessSession('basepath','D:\data\fNkx11\fNkx11_201110_sess19','analogChannelsList',65,'digitalChannelsList',[],'manualThr',true);
preprocessSession('basepath','D:\data\fNkx11\fNkx11_201103_sess14','analogChannelsList',65,'digitalChannelsList',[],'manualThr',true);


%% 03/18/22
updateExpFolder({'V:\data\fCck3', 'Y:\fCck3'},'Z:\data\fCck3');
batch_preprocessSession('basepath','Z:\data\fCck3','analysisPath','C:\data\fCck3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

updateExpFolder({'V:\data\fCr3', 'Y:\fCr3'},'Z:\data\fCr3');
batch_preprocessSession('basepath','Z:\data\fCr3','analysisPath','C:\data\fCr3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

updateExpFolder({'X:\fCck4', 'Y:\fCck4'},'Z:\data\fCck4');
batch_preprocessSession('basepath','Z:\data\fCck4','analysisPath','C:\data\fCck4','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
