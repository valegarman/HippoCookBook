%% WT5
updateExpFolder_temp({'D:\wt5'},'C:\wt5');

preprocessSession('basepath','C:\wt5\wt5_250317_sess11','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\wt5\wt5_250324_sess16','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);


processSession('basepath','C:\wt5\wt5_250314_sess10','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

%% WT6
processSession('basepath','Y:\wt6\wt6_250306_sess2','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','Y:\wt6\wt6_250307_sess3','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','Y:\wt6\wt6_250310_sess4','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','Y:\wt6\wt6_250306_sess2');
indexNewSession('basepath','Y:\wt6\wt6_250307_sess3');
indexNewSession('basepath','Y:\wt6\wt6_250310_sess4');



%% WT1
updateExpFolder_temp({'D:\wt1'},'C:\wt1');
preprocessSession('basepath','C:\wt1\wt1_241205_sess15','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
<<<<<<< HEAD
=======


>>>>>>> 2262c805a068d48770f291c0eb9783831492ee25
%% WT3

updateExpFolder_temp({'D:\wt3'},'C:\wt3');

preprocessSession('basepath','C:\wt3\wt3_241122_sess11','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\wt3\wt3_241125_sess12','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\wt3\wt3_241126_sess13','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\wt3\wt3_241205_sess20','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);


batch_preprocessSession('basepath','C:\wt3','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);


processSession('basepath','C:\wt3\wt3_241120_sess9','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\wt3\wt3_241125_sess12','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

%% fCamk8

preprocessSession('basepath','C:\fCamk8\fCamk8_220609_sess8','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',[],'getPos',false,'sessionSummary',false)
% batch_preprocessSession('basepath','C:\wt3','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% fAxo2
preprocessSession('basepath','Y:\others_databases\fAxo2_180412_sess1','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',[1],'getPos',false,'sessionSummary',false)
preprocessSession('basepath','Y:\others_databases\fAxo2_180414_sess2','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',[1],'getPos',false,'sessionSummary',false)
preprocessSession('basepath','Y:\others_databases\fAxo3_180415_sess1','analysisPath',[],'cleanArtifacts',({[],[1]}),'analogChannelsList',[],'digitalChannelsList',[1],'getPos',false,'sessionSummary',false)

processSession('basepath','C:\fAxo2\fAxo2_180414_sess2','digital_optogenetic_channels',[1],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\fAxo2\fAxo2_180414_sess2');

%% fAxo3
processSession('basepath','C:\fAxo3\fAxo3_180415_sess1','digital_optogenetic_channels',[1],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
indexNewSession('basepath','C:\fAxo3\fAxo3_180415_sess1');

