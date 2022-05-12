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




try
    updateExpFolder({'X:\fCr1', 'Y:\fCr1'},'Z:\data\fCr1');
    batch_preprocessSession('basepath','Z:\data\fCr1','analysisPath','C:\data\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
end
% updateExpFolder({'Y:\fCamk7','U:\fCamk7'},'J:\fCamk7');
% batch_preprocessSession('basepath','J:\fCamk7','analysisPath','C:\data\fCamk7','cleanArtifacts',({[3 4 5 6 7 8],[1 8 11 12 13 14 15 16 17]}),'analogChannelsList',[3 4 5 6 7 8],'digitalChannelsList',[1 8 11 12 13 14 15 16 17]);
updateExpFolder({'V:\data\fCr3', 'Y:\fCr3'},'Z:\data\fCr3');

updateExpFolder({'V:\data\fCr4','Y:\fCr4'},'Z:\data\fCr4');

batch_preprocessSession('basepath','Z:\data\fCr4','analysisPath','C:\data\fCr4','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

updateExpFolder('V:\data\fCr4','Z:\data\fCr4');
batch_preprocessSession('basepath','Z:\data\fCr4','analysisPath','C:\data\fCr4','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);
preprocessSession('basepath','C:\data\fCr4\fCr4_220418_sess9','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',true);




%% fPV3
batch_preprocessSession('basepath','J:\fPv3','analysisPath','C:\data\fPv3','cleanArtifacts',({1,[]}),'analogChannelsList',1,'digitalChannelsList',[],'bazler_ttl_channel',1);
preprocessSession('basepath','C:\data\fPv4\fPv4_210311_sess4','analogChannelsList',1,'digitalChannelsList',[],'bazler_ttl_channel',1);

try
    batch_sessionSummary('basepath','J:\fPv3','analogChannelsList',1,'digitalChannelsList',[]);
end

try
    batch_sessionSummary('basepath','J:\fPv4','analogChannelsList',1,'digitalChannelsList',[]);
end

batch_preprocessSession('basepath','J:\fPv4','analysisPath','C:\data\fPv4','cleanArtifacts',({1,[]}),'analogChannelsList',1,'digitalChannelsList',[],'bazler_ttl_channel',1);

%% fCamk7
preprocessSession('basepath','J:\fCamk7\fCamk7_220420_sess16','exclude','analogPulses','analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10);


batch_preprocessSession('basepath','J:\fCamk7','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',1);

preprocessSession('basepath','J:\fCamk7\fCamk7_220418_sess14','analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220419_sess15','analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);


preprocessSession('basepath','J:\fCamk7\fCamk7_220421_sess17','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220422_sess18','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);


computeSessionSummary('basepath','J:\fCamk7\fCamk7_220418_sess14','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220419_sess15','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220421_sess17','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220422_sess18','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

updateExpFolder({'Y:\fCamk7','T:\fCamk7','X:\fCamk7','T:\'},'J:\fCamk7');
preprocessSession('basepath','J:\fCamk7\fCamk7_220425_sess19','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220426_sess20','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220427_sess21','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220428_sess22','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220429_sess23','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);

preprocessSession('basepath','J:\fCamk7\fCamk7_220504_sess24','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220505_sess25','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220506_sess26','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220509_sess27','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220510_sess28','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220511_sess29','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);



computeSessionSummary('basepath','J:\fCamk7\fCamk7_220425_sess19','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220426_sess20','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220427_sess21','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

computeSessionSummary('basepath','J:\fCamk7\fCamk7_220429_sess23','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

updateExpFolder({'Y:\fCamk7','T:\fCamk7','X:\fCamk7'},'J:\fCamk7');





computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220302_sess23','analogChannelsList',[],'digitalChannelsList',[1]);
computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220303_sess24','analogChannelsList',[],'digitalChannelsList',[1]);
computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220304_sess25','analogChannelsList',[],'digitalChannelsList',[1]);



getuLEDsPulses();



updateExpFolder({'Y:\fCamk7','T:\fCamk7','X:\fCamk7'},'J:\fCamk7');
updateExpFolder({'Y:\fCamk7','T:\fCamk7'},'J:\fCamk7');

editDatFile(pwd,[0 839]);


updateExpFolder('V:\data\fId2','Z:\data\fId2');
batch_preprocessSession('basepath','Z:\data\fId2','analysisPath','C:\data\fId2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',1);
computeSessionSummary('basepath','Z:\data\fId2\fId2_220510_sess1','analogChannelsList',[],'digitalChannelsList',1);

