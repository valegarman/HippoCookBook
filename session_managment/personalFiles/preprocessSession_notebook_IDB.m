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



%% Once they have also been processed, then:
% 
% indexNewSession('basepath','Z:\astro7\astro7_260210_sess1');
% indexNewSession('basepath','Z:\astro7\astro7_260211_sess2');
% indexNewSession('basepath','Z:\astro7\astro7_260214_sess5');
% indexNewSession('basepath','Z:\astro7\astro7_260217_sess7');
% indexNewSession('basepath','Z:\astro7\astro7_260222_sess12');
% 
% indexNewSession('basepath','Z:\astro8\astro8_260212_sess1');
% indexNewSession('basepath','Z:\astro8\astro8_260213_sess2');
% indexNewSession('basepath','Z:\astro8\astro8_260214_sess3');
% indexNewSession('basepath','Z:\astro8\astro8_260216_sess4');
% indexNewSession('basepath','Z:\astro8\astro8_260219_sess7');
indexNewSession('basepath','Z:\astro8\astro8_260220_sess8');
indexNewSession('basepath','Z:\astro8\astro8_260221_sess9');
indexNewSession('basepath','Z:\astro8\astro8_260224_sess12');
indexNewSession('basepath','Z:\astro8\astro8_260228_sess16');
indexNewSession('basepath','Z:\astro8\astro8_260302_sess18');


