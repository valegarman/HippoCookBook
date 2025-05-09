%% APP2
updateExpFolder_temp({'E:\app2'},'C:\app2');

%preprocessSession('basepath','C:\app2\app2_250331_sess1','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

processSession('basepath','C:\app2\app2_250401_sess2','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\wt5\wt5_250320_sess14');

%% WT5
preprocessSession('basepath','C:\wt5\wt5_250331_sess1','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
