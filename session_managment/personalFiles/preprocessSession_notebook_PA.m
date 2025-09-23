%% KOBE
preprocessSession('basepath','D:\kobe\kobe_000001_sess1','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);



%% WT9
updateExpFolder_temp({'E:\wt9'},'D:\wt9');

preprocessSession('basepath','D:\wt9\wt9_250917_sess3','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\wt9\wt9_250917_sess3Opto','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% CANCER4
updateExpFolder_temp({'E:\cancer4'},'D:\cancer4');

preprocessSession('basepath','D:\cancer4\cancer4_250722_sess7','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\cancer4\cancer4_250724_sess9','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% CANCER5
updateExpFolder_temp({'E:\cancer5'},'D:\cancer5');

preprocessSession('basepath','D:\cancer5\cancer5_250719_sess4','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

processSession('basepath','D:\cancer5\cancer5_250719_sess4','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\cancer5\cancer5_250721_sess6','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\data\cancer5\cancer5_250716_sess1');
indexNewSession('basepath','C:\data\cancer5\cancer5_250721_sess6');

%% CANCER6
updateExpFolder_temp({'E:\cancer6'},'D:\cancer6');

%% ASTRO5
updateExpFolder_temp({'E:\astro5'},'F:\astro5');

preprocessSession('basepath','F:\astro5\astro5_250728_sess20','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
%% ASTRO6
updateExpFolder_temp({'E:\astro6'},'F:\astro6');

preprocessSession('basepath','F:\astro6\astro6_250728_sess20','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% CANCER2
updateExpFolder_temp({'E:\cancer2'},'D:\cancer2');
batch_preprocessSession('basepath','D:\cancer2','sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\cancer2\cancer2_250708_sess7','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

processSession('basepath','Y:\cancer2\cancer2_250710_sess9','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\data\cancer2\cancer2_250708_sess7');
indexNewSession('basepath','Y:\cancer2\cancer2_250710_sess9');

%% CANCER 3
updateExpFolder_temp({'E:\cancer3'},'D:\cancer3');

preprocessSession('basepath','D:\cancer3\cancer3_250708_sess7','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
processSession('basepath','C:\data\cancer3\cancer3_250709_sess8','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\data\cancer3\cancer3_250703_sess2');
indexNewSession('basepath','C:\data\cancer3\cancer3_250708_sess7');
indexNewSession('basepath','C:\data\cancer3\cancer3_250709_sess8');

%% CAMK13
updateExpFolder_temp({'E:\camk13'},'D:\camk13');
preprocessSession('basepath','D:\camk13\camk13_250618_sess5','medianSubstr',[21 28 26 22 23 27 17 32 30 18 16 1 15 3 2 14 12 5 11],'analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[1 2]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

processSession('basepath','Y:\camk13\camk13_250618_sess5','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','Y:\camk13\camk13_250702_sess13','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250717_sess28','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250723_sess32','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','Y:\camk13\camk13_250722_sess31','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

processSession('basepath','Y:\camk13\camk13_250718_sess29','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250731_sess38','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250714_sess25','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250715_sess26','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250713_sess24','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\camk13\camk13_250712_sess23','digital_optogenetic_channels',[1 2],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','Y:\camk13\camk13_250618_sess5');
indexNewSession('basepath','Y:\camk13\camk13_250702_sess13');
indexNewSession('basepath','C:\data\camk13\camk13_250717_sess28');
indexNewSession('basepath','C:\data\camk13\camk13_250723_sess32');
indexNewSession('basepath','C:\data\camk13\camk13_250731_sess38');
indexNewSession('basepath','C:\data\camk13\camk13_250714_sess25');
indexNewSession('basepath','C:\data\camk13\camk13_250715_sess26');
indexNewSession('basepath','C:\data\camk13\camk13_250713_sess24');
indexNewSession('basepath','C:\data\camk13\camk13_250712_sess23');

%% ASTRO3
updateExpFolder_temp({'E:\astro3'},'D:\astro3');
batch_preprocessSession('basepath','D:\astro3','sessionSummary',false,'getPos',false);

preprocessSession('basepath','D:\astro3\astro3_250529_sess1','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\astro3\astro3_250530_sess2','analysisPath',[],'exclude_shanks',3,'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\astro3\astro3_250604_sess5','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','D:\astro3\astro3_250617_sess18','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% APP2
updateExpFolder_temp({'D:\app2'},'E:\app2');

preprocessSession('basepath','C:\app2\app2_250331_sess1','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','E:\app2\app2_250415_sess12','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','E:\app2\app2_250416_sess13','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','F:\app2\app2_250411_sess10','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

preprocessSession('basepath','F:\app2\app2_250507_sess17','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','F:\app2\app2_250508_sess18','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

batch_preprocessSession('basepath','E:\app2','sessionSummary',false,'getPos',false);

processSession('basepath','C:\app2\app2_250331_sess1','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\app2\app2_250411_sess10','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\app2\app2_250415_sess12','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\data\app2\app2_250411_sess10');



%% WT7
updateExpFolder_temp({'E:\wt7'},'D:\wt7');

preprocessSession('basepath','E:\wt7\wt7_250415_090618','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
batch_preprocessSession('basepath','D:\wt7','sessionSummary',false,'getPos',false);

editDatFile('E:\wt7\wt7_250404_sess1\wt7_250404_144632',[972 1018],'option','remove');

batch_preprocessSession('basepath','E:\wt7','sessionSummary',false,'getPos',false);

processSession('basepath','C:\data\wt7\wt7_250404_sess1','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','C:\data\wt7\wt7_250407_sess2','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','E:\wt7\wt7_250408_sess3','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','E:\wt7\wt7_250409_sess4','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','D:\wt7 sorted\wt7_250414_sess7','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);


indexNewSession('basepath','C:\data\wt7\wt7_250404_sess1');
indexNewSession('basepath','C:\data\wt7\wt7_250407_sess2');
indexNewSession('basepath','E:\wt7\wt7_250408_sess3');
indexNewSession('basepath','E:\wt7\wt7_250409_sess4');
indexNewSession('basepath','D:\wt7 sorted\wt7_250414_sess7');

%% WT5
updateExpFolder_temp({'D:\wt5'},'C:\wt5');

preprocessSession('basepath','C:\wt5\wt5_250317_sess11','analysisPath',[],'exclude_shanks',[3],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

processSession('basepath','C:\wt5\wt5_250326_sess18','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','C:\wt5\wt5_250326_sess18');
indexNewSession('basepath','Y:\wt5\wt5_250320_sess14');
%% WT6
processSession('basepath','Y:\wt6\wt6_250306_sess2','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

indexNewSession('basepath','Y:\wt6\wt6_250310_sess4');

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


processSession('basepath','Y:\wt3\wt3_241119_sess8','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);
processSession('basepath','Y:\wt3\wt3_241120_sess9','digital_optogenetic_channels',[],'analog_optogenetic_channels',[],'promt_hippo_layers',true);

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

