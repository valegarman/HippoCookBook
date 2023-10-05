preprocessSession('basepath','Z:\data\fNkx11\fNkx11_201030_sess12b','analysisPath','C:\data\fNkx11','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',[],'bazler_ttl_channel',1);

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
computeSessionSummary('basepath','W:\Buzsakilabspace\Datasets\ValeroM\fCck3\fCck3_220308_sess13','analogChannelsList',[],'digitalChannelsList',1);


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

% fCr2
batch_preprocessSession('basepath','Z:\data\fCr2','analysisPath','C:\data\fCr2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','Z:\data\fCr2','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
computeSessionSummary('basepath','Z:\data\fCr2\fCr2_220306_sess5','analogChannelsList',[],'digitalChannelsList',1);



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
preprocessSession('basepath','C:\data\fCr4\fCr4_220429_sess17','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',true);




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
preprocessSession('basepath','J:\fCamk7\fCamk7_220425_sess19','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);

preprocessSession('basepath','Z:\data\fCamk7\fCamk7_220425_sess19','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);

computeSessionSummary('basepath','J:\fCamk7\fCamk7_220418_sess14','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220419_sess15','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220421_sess17','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220422_sess18','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);


updateExpFolder({'Y:\fCamk7','T:\fCamk7'},'J:\fCamk7');

preprocessSession('basepath','J:\fCamk7\fCamk7_220425_sess19','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220426_sess20','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220427_sess21','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);



preprocessSession('basepath','J:\fCamk7\fCamk7_220505_sess25','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220506_sess26','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220509_sess27','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220510_sess28','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220511_sess29','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);


preprocessSession('basepath','J:\fCamk7\fCamk7_220425_sess19','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220426_sess20','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220427_sess21','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220428_sess22','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);
preprocessSession('basepath','J:\fCamk7\fCamk7_220429_sess23','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);

preprocessSession('basepath','J:\fCamk7\fCamk7_220504_sess24','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220505_sess25','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220506_sess26','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220509_sess27','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false);

preprocessSession('basepath','J:\fCamk7\fCamk7_220510_sess28','cleanArtifacts',[],'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220511_sess29','cleanArtifacts',[],'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false);
preprocessSession('basepath','J:\fCamk7\fCamk7_220514_sess30','cleanArtifacts',[],'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false);

preprocessSession('basepath','D:\fCamk7\fCamk7_220514_sess30','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',true);


computeSessionSummary('basepath','J:\fCamk7\fCamk7_220425_sess19','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220426_sess20','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
computeSessionSummary('basepath','J:\fCamk7\fCamk7_220427_sess21','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

computeSessionSummary('basepath','J:\fCamk7\fCamk7_220504_sess24','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

computeSessionSummary('basepath','J:\fCamk7\fCamk7_220429_sess23','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);

updateExpFolder({'Y:\fCamk7','T:\fCamk7','X:\fCamk7'},'J:\fCamk7');





computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220302_sess23','analogChannelsList',[],'digitalChannelsList',[1]);
computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220303_sess24','analogChannelsList',[],'digitalChannelsList',[1]);
computeSessionSummary('basepath','Z:\data\fCck1\fCck1_220304_sess25','analogChannelsList',[],'digitalChannelsList',[1]);

computeSessionSummary('basepath','W:\Buzsakilabspace\Datasets\ValeroM\fCck1\fCck1_220131_sess1','analogChannelsList',1,'digitalChannelsList',[]);


getuLEDsPulses();



updateExpFolder({'Y:\fCamk7','T:\fCamk7','X:\fCamk7'},'J:\fCamk7');
updateExpFolder({'Y:\fCamk7','T:\fCamk7'},'J:\fCamk7');

editDatFile(pwd,[0 839]);
computeSessionSummary('basepath','K:\fCamk7\fCamk7_220511_sess29','analogChannelsList',[],'digitalChannelsList',[1]);

%% fId2
updateExpFolder('V:\data\fId2','X:\data\fId2');
% batch_preprocessSession('basepath','Z:\data\fId2','analysisPath','C:\data\fId2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','Z:\data\fId2','analysisPath','C:\data\fId2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fVIP1
updateExpFolder({'V:\data\fVIP1','Y:\fVIP1'},'Z:\data\fVIP1');
% batch_preprocessSession('basepath','Z:\data\fVIP1','analysisPath','C:\data\fVIP1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','Z:\data\fVIP1','analysisPath','C:\data\fVIP1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fId3
updateExpFolder('V:\data\fId3','Z:\data\fId3');
% batch_preprocessSession('basepath','Z:\data\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','Z:\data\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fNkx12
updateExpFolder('V:\data\fNkx12','X:\data\fNxk12');
% batch_preprocessSession('basepath','Z:\data\fNxk12','analysisPath','C:\data\fNxk12','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','X:\data\fNxk12','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);


%% fCr2
batch_preprocessSession('basepath','X:\data\fCr2','analysisPath','C:\data\fCr2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
batch_preprocessSession('basepath','X:\data\fCr2','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);
% updateExpFolder('V:\data\fId2','Z:\data\fId2');
% batch_preprocessSession('basepath','Z:\data\fId2','analysisPath','C:\data\fId2','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fVIP1
updateExpFolder({'V:\data\fVIP1','Y:\fVIP1'},'Z:\data\fVIP1');
batch_preprocessSession('basepath','Z:\data\fVIP1','analysisPath','C:\data\fVIP1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fId3
updateExpFolder('V:\data\fId3','Z:\data\fId3');
batch_preprocessSession('basepath','Z:\data\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fNkx12
updateExpFolder('T:\fNkx12','Z:\data\fNkx12');
batch_preprocessSession('basepath','Z:\data\fNkx12','analysisPath','C:\data\fNkx12','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fCck1
batch_preprocessSession('basepath','W:\Buzsakilabspace\Datasets\ValeroM\fCck1','analysisPath','C:\data\fCck1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fCr3
batch_preprocessSession('basepath','Z:\data\fCr3','analysisPath','C:\data\fCr3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fCck4
batch_preprocessSession('basepath','J:\fCck4','analysisPath','C:\data\fCck4','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);


basepath = 'Z:\data\fVIP1\fVIP1_220526_sess12';
cd(basepath);
optogeneticResponses = getOptogeneticResponse;



%%
preprocessSession('basepath','J:\fCamk7\fCamk7_220509_sess27','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

%%
batch_preprocessSession('basepath','Z:\data\fCr3','analysisPath','C:\data\fCr3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'getPos',false);
batch_preprocessSession('basepath','W:\Buzsakilabspace\Datasets\ValeroM\fCck3','analysisPath','C:\data\fCck3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'getPos',false);
batch_preprocessSession('basepath','W:\Buzsakilabspace\Datasets\ValeroM\fCck4','analysisPath','C:\data\fCck4','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'getPos',false);
%% fVIP1
updateExpFolder({'V:\data\fVIP1','Y:\fVIP1'},'Z:\data\fVIP1');
batch_preprocessSession('basepath','Z:\data\fVIP1','analysisPath','C:\data\fVIP1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fId3
updateExpFolder('V:\data\fId3','Z:\data\fId3');
batch_preprocessSession('basepath','Z:\data\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);

%% fPv5
updateExpFolder({'T:\fPv5','V:\data\fPv5'},'Z:\data\fId3');
batch_preprocessSession('basepath','Z:\data\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);


%% fCamk8
updateExpFolder({'T:\fCamk8','Y:\fCamk8'},'J:\fCamk8');
batch_preprocessSession('basepath','J:\fCamk8','cleanArtifacts',[],'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false);

%%

updateExpFolder({'V:\data\fVIP1','Y:\fVIP1'},'Z:\data\fVIP1');
updateExpFolder({'T:\fPv5','V:\data\fPv5'},'Z:\data\fPv5');
updateExpFolder({'T:\fCamk8','Y:\fCamk8'},'J:\fCamk8');

batch_preprocessSession('basepath','Z:\data\fVIP1','analysisPath','C:\data\fVIP1','cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',10,'getPos',false);



%% FL4
changeFilesName('Z:\FLR\FL4\FL4_080322');
changeFilesName('Z:\FLR\FL4\FL4_090322');
changeFilesName('Z:\FLR\FL4\FL4_100322');
changeFilesName('Z:\FLR\FL4\FL4_110322','socialParadigm',true);
changeFilesName('basepath','Z:\FLR\FL4\FL4_170322','generalPath','Z:\FLR','socialParadigm',true);
changeFilesName('basepath','Z:\FLR\FL4\FL4_180322','generalPath','Z:\FLR','socialParadigm',true);
changeFilesName('basepath','Z:\FLR\FL4\FL4_210322','generalPath','Z:\FLR','socialParadigm',true);
changeFilesName('basepath','Z:\FLR\FL4\FL4_220322','generalPath','Z:\FLR','socialParadigm',true);
% updateExpFolder('Z:\FLR\FL4','D:\FLR\FL4');
basepath = 'D:\FLR\FL4';
cd(basepath);
arrangeSessionFolder;
createFiles('basepath',basepath);

basepath = 'D:\FLR\FL4\FL4_080322_sess1';
preprocessSession_pablo('basepath',basepath,'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true,'anyMaze',true);
digitalChannelsList = [];
analogChannelsList = [];
computeSessionSummary_pablo('digitalChannelsList',digitalChannelsList,'analogChannelsList',analogChannelsList);

batch_preprocessSession_pablo('basepath','D:\FLR\FL4','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

basepath = 'D:\FLR\FL4\FL4_080322_sess1';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);

basepath = 'D:\FLR\FL4\FL4_090322_sess2';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);

batch_sessionSummary_pablo('basepath','D:\FLR\FL4','analogChannelsList',[],'digitalChannelsList',[]);


preprocessSession_pablo('basepath','D:\FLR\FL3\FL3_020322_sess1','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_080322_sess1','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_090322_sess2','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_100322_sess3','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_110322_sess4','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_170322_sess5','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_180322_sess6','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_210322_sess7','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\FLR\FL4\FL4_220322_sess8','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);


%% FL3
batch_changeFilesName('basepath','Z:\FLR\FL3','generalPath','Z:\FLR','socialParadigm',true);
% updateExpFolder('Z:\FLR\FL5','D:\FLR\FL5');
basepath = 'D:\FLR\FL3';
cd(basepath);
arrangeSessionFolder;
createFiles('basepath',basepath);
preprocessSession_pablo('basepath','D:\FLR\FL3\FL3_020322_sess1','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
batch_preprocessSession_pablo('basepath','D:\FLR\FL3','analysisPath','C:\FL3','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

preprocessSession_pablo('basepath','D:\FLR\FL3\FL3_170322_sess6','analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

%% FL5
batch_changeFilesName('basepath','Z:\FLR\FL5','generalPath','Z:\FLR','socialParadigm',true);
% updateExpFolder('Z:\FLR\FL3','D:\FLR\FL3');
basepath = 'D:\FLR\FL5';
cd(basepath);
arrangeSessionFolder;
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\FLR\FL5','analysisPath','C:\FL5','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);


%% FL9
basepath = 'D:\FLR\FL9';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\FLR\FL5','analysisPath','C:\FL5','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);


%% HPS22
basepath = 'F:\data\HPS22';
basepath = 'D:\HPS22';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\HPS22','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
computeSessionSummary_pablo('basepath','F:\data\HPS22\HPS22_040621_sess25','analogChannelsList',[],'digitalChannelsList',[]);
computeSessionSummary_pablo('basepath','F:\data\HPS22\HPS22_140521_sess13','analogChannelsList',[],'digitalChannelsList',[]);

computeSessionSummary_pablo('basepath','F:\data\HPS22\HPS22_040621_sess25','analogChannelsList',[],'digitalChannelsList',[]);
batch_sessionSummary_pablo('basepath','F:\data\HPS22','analogChannelsList',[],'digitalChannelsList',[]);


preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_140521_sess13','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_270421_sess7','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_230421_sess6','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_290421_sess9','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_130521_sess11','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_270521_sess18','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_280521_sess19','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_040621_sess125','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)

preprocessSession_pablo('basepath','D:\HPS22\HPS22_130521_sess11','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_040621_sess125','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)

preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_210421_sess4','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false)
preprocessSession_pablo('basepath','F:\data\HPS22\HPS22_220421_sess5','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false)

preprocessSession_pablo('basepath','D:\HPS22\HPS22_230421_sess6','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false)


batch_preprocessSession_pablo('basepath','D:\HPS22','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

%% HPS23
basepath = 'F:\data\HPS23';
basepath = 'D:\HPS23';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\HPS23','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

basepath = 'F:\data\HPS23\HPS23_090621_sess9';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);
preprocessSession_pablo('basepath','D:\HPS23\HPS23_190521_sess5','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','D:\HPS23\HPS23_200521_sess6','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
preprocessSession_pablo('basepath','D:\HPS23\HPS23_210521_sess7','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\HPS23\HPS23_270521_sess8','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','F:\data\HPS23\HPS23_110621_sess10','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

batch_preprocessSession_pablo('basepath','D:\HPS23','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
preprocessSession_pablo('basepath','D:\HPS23\HPS23_180521_sess16','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

%% HPS24
basepath = 'F:\data\HPS24';
basepath = 'D:\HPS24';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','F:\data\HPS24','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

preprocessSession_pablo('basepath','F:\data\HPS24\HPS24_280621_sess5','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
preprocessSession_pablo('basepath','F:\data\HPS24\HPS24_290621_sess6','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

batch_preprocessSession_pablo('basepath','D:\HPS24','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

%% HPS25
basepath = 'F:\data\HPS25';
basepath = 'D:\HPS25';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','F:\data\HPS25','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
preprocessSession_pablo('basepath','F:\data\HPS25\HPS25_280621_sess3','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',false,'sessionSummary',false);

batch_preprocessSession_pablo('basepath','D:\HPS25','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

%% IPO135
basepath = 'F:\data\IPO135';
createFiles('basepath',basepath);

%% IPO150
preprocessSession_pablo('basepath','F:\data\IPO150\IPO150_261021_sess2','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true);
preprocessSession_pablo('basepath','F:\data\IPO150\IPO150_291021_sess3','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true);

%% IPO4311
preprocessSession_pablo('basepath','F:\data\IPO431\IPO431_221021_sess1','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true);


%% HPR21409

basepath = 'I:\HPR21409\HPR21409_150921_sess1';
createFiles('basepath',basepath);



%%
computeSessionSummary_pablo('basepath','J:\buzsakilab data\fCr1_220402_sess26','analogChannelsList',[],'digitalChannelsList',[1]);
computeSessionSummary_pablo('basepath','J:\buzsakilab data\fCr1_220404_sess27','analogChannelsList',[],'digitalChannelsList',[1]);



%% HPR21409
basepath = 'I:\HPR21409';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','I:\HPR21409','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);
basepath = 'I:\HPR21409\HPR21409_230921_sess3';
preprocessSession_pablo('basepath',basepath,'cleanArtifacts',({[],[7]}),'analogChannelsList',[],'digitalChannelsList',[7],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true,'anyMaze',true);

batch_preprocessSession_pablo('basepath','I:\HPR21409','analysisPath',[],'cleanArtifacts',({[],[7]}),'analogChannelsList',[],'digitalChannelsList',[7],'anymaze_ttl_channel',2);

basepath = 'I:\HPR21409\HPR21409_051021_sess10';
preprocessSession_pablo('basepath',basepath,'cleanArtifacts',({[],[7]}),'analogChannelsList',[],'digitalChannelsList',[7],'anymaze_ttl_channel',2,'getPos',true,'sessionSummary',true,'anyMaze',true);


%% VB20R3
basepath = 'I:\VB20R3';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','I:\HPR21409','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);


%% fCamk7

basepath = 'K:\fCamk7\fCamk7_220511_sess29';
computeSessionSummary('basepath',basepath,'exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[]);

computeSessionSummary('basepath','K:\fCamk7\fCamk7_220510_sess28','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[]);
computeSessionSummary('basepath','K:\fCamk7\fCamk7_220514_sess30','exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[]);

%%
createProbe('excel_file','electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-4');
createProbe('excel_file','electrodes_coordinates_DiagnosticBiochip-128-6-128ch&uLED_12LED-32Ch-4Shanks');


%% HM
createProbe('excel_file','electrodes_coordinates_UtahArray-96ch');
selectProbe();
basepath = 'D:\HM';
cd(basepath);
basepath = 'D:\HM2';
createNSFiles('basepath',basepath,'nChannels',96);

preprocessSession_pablo('basepath','D:\HM\HM_290622_sess1','analysisPath',[],'cleanArtifacts',({[1 2],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\HM\HM_291122_sess3','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[1],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\HM\HM_391122_sess4','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[1],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);

computeSessionSummary_pablo('basepath','F:\data\HM\HM_270722_sess2','analogChannelsList',[1 2 3],'digitalChannelsList',[]);


preprocessSession('basepath','K:\fCamk7\fCamk7_220511_sess29','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession('basepath','K:\fCamk7\fCamk7_220511_sess29','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);


createFiles('basepath','D:\FLR\FL10');
preprocessSession_pablo('basepath','D:\HM2\HM2_210223_sess1','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);



%% fSUB1
preprocessSession_pablo('basepath','J:\fSUB1\fSUB1_181122_sess7','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','D:\fSUB1\fSUB1_150223_sess37','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_110223_sess33','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_060223_sess28','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_120223_sess34','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_100223_sess32','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_261222_sess25','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','C:\DATA\fSUB1\fSUB1_291122_sess14','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);


%% fSUB2
preprocessSession_pablo('basepath','J:\fSUB2\fSUB2_150223_sess6','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','J:\fSUB2\fSUB2_160223_sess7','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);


%% HM2
bpath = 'C:\DATA\HM2';
createNSFiles('basepath',bpath,'nChannels',96);
preprocessSession_pablo('basepath','C:\DATA\HM2\HM2_090223_sess3','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','J:\fSUB2\fSUB2_220223_sess10','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);


%% HMX
basepath = 'J:\HM2';
createNSFiles('basepath',basepath,'nChannels',96);
cd(basepath);
preprocessSession_pablo('basepath','J:\HM2\HM2_270323_sess2','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);
preprocessSession_pablo('basepath','J:\HM2\HM2_160323_sess1','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','C:\DATA\HM2\HM2_28042_sess4','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);



%%

preprocessSession_pablo('basepath','C:\DATA\IPO11700_280723_sess19','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

% preprocessSession_pablo('basepath','J:\fSUB1\fSUB1_181122_sess7','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'changeAnalogInputs',false,'getPos',false,'sessionSummary',false);
% 
%% C:\DATA\IPO14369_260723_sess17 

preprocessSession_pablo('basepath','D:\IPO14369_270623_sess18','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

preprocessSession_pablo('basepath','D:\Recordings\IPO14369\IPO14369_030823_sess20','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'bazler_ttl_channel',10,'getPos',false,'sessionSummary',false);

