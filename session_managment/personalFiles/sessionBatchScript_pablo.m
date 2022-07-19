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
batch_preprocessSession_pablo('basepath','D:\FLR\FL4','analysisPath','C:\FL4','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

basepath = 'D:\FLR\FL4\FL4_080322_sess1';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);

basepath = 'D:\FLR\FL4\FL4_090322_sess2';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);

batch_sessionSummary_pablo('basepath','D:\FLR\FL4','analogChannelsList',[],'digitalChannelsList',[]);

%% FL3
batch_changeFilesName('basepath','Z:\FLR\FL3','generalPath','Z:\FLR','socialParadigm',true);
% updateExpFolder('Z:\FLR\FL5','D:\FLR\FL5');
basepath = 'D:\FLR\FL3';
cd(basepath);
arrangeSessionFolder;
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\FLR\FL3','analysisPath','C:\FL3','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

%% FL5
batch_changeFilesName('basepath','Z:\FLR\FL5','generalPath','Z:\FLR','socialParadigm',true);
% updateExpFolder('Z:\FLR\FL3','D:\FLR\FL3');
basepath = 'D:\FLR\FL5';
cd(basepath);
arrangeSessionFolder;
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','D:\FLR\FL5','analysisPath','C:\FL5','cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);


%% HPS22

%% HPS23
basepath = 'F:\data\HPS23';
createFiles('basepath',basepath);
batch_preprocessSession_pablo('basepath','F:\data\HPS23','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2);

basepath = 'F:\data\HPS23\HPS23_090621_sess9';
computeSessionSummary_pablo('basepath',basepath,'analogChannelsList',[],'digitalChannelsList',[]);

%% HPR21409

basepath = 'I:\HPR21409\HPR21409_150921_sess1';
createFiles('basepath',basepath);







