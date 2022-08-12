%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'U:\fCr1', 'Y:\fCr1'},'X:\data\fCr1');

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

updateExpFolder('V:\data\fId3','J:\fId3');


%2% Preprocessing
batch_preprocessSession('basepath','X:\data\fCr1','analysisPath','F:\fCr1','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1);

batch_preprocessSession('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

batch_preprocessSession('basepath','D:\Dropbox\DATA\anna_data\NY14');

batch_preprocessSession('basepath','J:\fId3','analysisPath','C:\data\fId3','cleanArtifacts',({[],1}),'analogChannelsList',[],'digitalChannelsList',1,'bazler_ttl_channel',[]);

preprocessSession('basepath',pwd,'analogChannelsList',analogChannelsList,'spikeSort',spikeSort,'getPos',getPos, 'cleanArtifacts',cleanArtifacts,...
                    'medianSubstr',medianSubstr,'tracking_pixel_cm',tracking_pixel_cm,'sessionSummary',sessionSummary,'digitalChannelsList',digitalChannelsList,'bazler_ttl_channel',bazler_ttl_channel);


=======
batch_preprocessSession('basepath','E:\data\NY13')

%3% Computing summary
batch_sessionSummary('basepath','J:\fId3','analogChannelsList',[],'digitalChannelsList',1);

batch_sessionSummary('basepath','G:\data\fNkx11','analogChannelsList',65,'digitalChannelsList',[]);


% After Spike sorting, got
processSession('analog_optogenetic_channels',65,'promt_hippo_layers',true,'manual_analog_pulses_threshold',true,'bazler_ttl_channel',1);
 
%% others
editDatFile(pwd,[26460 inf],'option','zeroes');   % remove zeroes
