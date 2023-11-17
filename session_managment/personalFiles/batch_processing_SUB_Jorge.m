


bpath = 'C:\DATA\fSUB1\fSUB1_130223_sess35'
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',21,'SWChannel',13,'thetaChannel',21,'excludeAnalysis',{'getHippocampalLayers'});


bpath = 'C:\DATA\fSUB1\fSUB1_140223_sess36'
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',22,'SWChannel',13,'thetaChannel',22,'excludeAnalysis',{'getHippocampalLayers'});
bpath = 'C:\DATA\fSUB1\fSUB1_150223_sess37'
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',22,'SWChannel',13,'thetaChannel',22,'excludeAnalysis',{'getHippocampalLayers'});

createNSFiles('basepath','C:\DATA\HM2')
bpath='C:\DATA\HM2\HM2_190523_sess5'
% preprocessSession_pablo('basepath','F:\data\IPO149\IPO149_281021_sess5','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
processSession_HM('basepath',bpath,'project','HMProject');
 batch_preprocessSession_pablo('basepath',bpath,'project','HMProject')
batch_preprocessSession('basepath','C:\DATA\HM2\HM2_190523_sess5','analysisPath',[])%,[],'cleanArtifacts',({[],[]})),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)


% createNSFiles('basepath','D:\HM');
% bpath = 'F:\data\HM\HM_270722_sess2'; 
% processSession_HM('basepath',bpath,'project','HMProject');
% indexNewSession_pablo('basepath',bpath,'project','HMProject');