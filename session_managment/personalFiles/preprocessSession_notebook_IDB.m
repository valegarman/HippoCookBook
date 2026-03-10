%% TRAUMA1
updateExpFolder_temp({'D:\trauma1'},'F:\trauma1');

preprocessSession('basepath','F:\trauma1\trauma1_260205_sess22','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','F:\trauma1\trauma1_260206_sess23','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','F:\trauma1\trauma1_260209_sess24','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);



%% Astro7
updateExpFolder_temp({'E:\astro77'},'C:\data\astro7');

% batch_preprocessSession('basepath','F:\astro7','cleanArtifacts',({[],[8]}),'sessionSummary',true);

% preprocessSession('basepath','C:\data\astro7\astro7_260212_sess3','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\data\astro7\astro7_260302_sess18','analysisPath',[],'exclude_shanks',[1],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
% preprocessSession('basepath','F:\astro7\astro7_260222_sess12','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);

%% Astro8
updateExpFolder_temp({'E:\astro87'},'C:\data\astro8');

% editDatFile(pwd,[17317, 17340],'option', 'remove')

% batch_preprocessSession('basepath','F:\astro7','cleanArtifacts',({[],[8]}),'sessionSummary',true);

% preprocessSession('basepath','C:\data\astro8\astro8_260212_sess1','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
preprocessSession('basepath','C:\data\astro8\astro8_260302_sess18','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
% preprocessSession('basepath','F:\astro8\astro8_260222_sess10','analysisPath',[],'exclude_shanks',[],'cleanArtifacts',({[],[8]}),'analogChannelsList',[],'digitalChannelsList',[],'sessionSummary',false,'getPos',false);
