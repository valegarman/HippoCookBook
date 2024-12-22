
%% WT1
updateExpFolder_temp({'D:\wt1'},'C:\wt1');
preprocessSession('basepath','C:\wt1\wt1_241205_sess15','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);



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

