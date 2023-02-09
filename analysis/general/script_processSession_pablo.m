%% SOCIAL PROJECT
createProbe('excel_file','electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-4');
createProbe('excel_file','electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-3');


% FL4_080322_sess1
bpath = 'D:\FLR\FL4\FL4_080322_sess1';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_090322 sess2
bpath = 'D:\FLR\FL4\FL4_090322_sess2';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_100322_sess3
bpath = 'D:\FLR\FL4\FL4_100322_sess3';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_110322_sess4
bpath = 'D:\FLR\FL4\FL4_110322_sess4';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_170322_sess5
bpath = 'D:\FLR\FL4\FL4_170322_sess5';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_180322_sess6
bpath = 'D:\FLR\FL4\FL4_180322_sess6';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_210322_sess7
bpath = 'D:\FLR\FL4\FL4_210322_sess7';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL4_220322_sess8
bpath = 'D:\FLR\FL4\FL4_220322_sess8';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_020322_sess1
bpath = 'D:\FLR\FL3\FL3_020322_sess1';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_030322_sess2
bpath = 'D:\FLR\FL3\FL3_030322_sess2';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_100322_sess3
bpath = 'D:\FLR\FL3\FL3_100322_sess3';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_110322_sess4
bpath = 'D:\FLR\FL3\FL3_110322_sess4';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_160322_sess5
bpath = 'D:\FLR\FL3\FL3_160322_sess5';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');

% FL3_170322_sess6
bpath = 'D:\FLR\FL3\FL3_170322_sess6';
processSession_pablo('basepath',bpath,'project','SocialProject','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',28,'SWChannel',4,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');
indexNewSession('basepath',bpath,'project','SocialProject');


%% SPATIAL SUBICULUM PROJECT
createProbe('excel_file','electrodes_coordinates_Buzsaki64(64 ch, 8 shanks, staggered)');
createProbe('excel_file','electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177-Allego');

% fSUB1_011122_sess16
bpath = 'D:\fSUB1\fSUB1_011222_sess16';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2);
indexNewSession_pablo('basepath',bpath,'project','SUBProject');

% fSUB1_121222_sess18 (Open Field)
bpath = 'D:\fSUB1\fSUB1_121222_sess18';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',47,'SWChannel',45,'thetaChannel',47);
indexNewSession_pablo('basepath',bpath,'project','SUBProject');

% fSUB1_231222_sess25 (TMaze)
bpath = 'D:\fSUB1\fSUB1_231222_sess25';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',18,'SWChannel',45,'thetaChannel',18);
indexNewSession_pablo('basepath',bpath,'project','SUBProject');


% HPS22_210521_sess17
bpath = 'F:\data\HPS22\HPS22_210521_sess17';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',39,'SWChannel',13,'thetaChannel',39,'excludeAnalysis',{'getHippocampalLayers'});
indexNewSession_pablo('basepath',bpath,'project','Subiculum Project');

% HPS22_210521_sess117 (to be done)
bpath = 'F:\data\HPS22\HPS22_210521_sess117';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',39,'SWChannel',13,'thetaChannel',39,'excludeAnalysis',{'getHippocampalLayers'});
indexNewSession_pablo('basepath',bpath,'project','SubiculumProject');



% fSUB1_281122_sess13 (Open Field)
bpath = 'D:\fSUB1\fSUB1_281122_sess13';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'tint',true,'speedThresh',0,'gridAnalysis',true,'randomization',true);
indexNewSession_pablo('basepath',bpath,'project','SubiculumProject');




%% MK801 Project

% HPS22_100621_sess26 MK801
bpath = 'F:\data\HPS22\HPS22_100621_sess26';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',24,'SWChannel',13,'thetaChannel',24);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS22_150621_sess27 VEHICLE
bpath = 'F:\data\HPS22\HPS22_150621_sess27';
processSession_pablo('basepath',bpath,'promt_hippo_layers',false,'anymaze_ttl_channel',2,'rippleChannel',24,'SWChannel',13,'thetaChannel',24);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS22_160621_sess28 KETAMINE (Animal dead when injected ketamine. Not run
% processSession)
bpath = 'F:\data\HPS22\HPS22_160621_sess28';
processSession_pablo('basepath',bpath,'promt_hippo_layers',false,'anymaze_ttl_channel',2,'rippleChannel',24,'SWChannel',13,'thetaChannel',24);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS23_090621_sess9 MK801
bpath = 'F:\data\HPS23\HPS23_090621_sess9';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',30,'SWChannel',13,'thetaChannel',30);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS23_110621_sess210 VEHICLE 
bpath = 'F:\data\HPS23\HPS23_110621_sess10'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',30,'SWChannel',13,'thetaChannel',30);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS23_160621_sess11 KETAMINE
bpath = 'F:\data\HPS23\HPS23_160621_sess11';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',30,'SWChannel',13,'thetaChannel',30);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS24_280621_sess5 VEHICLE ( I don't understand the recording. Probably
% headstage connected on the contrary ??)
bpath = 'F:\data\HPS24\HPS24_280621_sess5';
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS24_290621_sess6 MK801
bpath = 'F:\data\HPS24\HPS24_290621_sess6';
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS24_230621_sess4 KETAMINE 
bpath = 'F:\data\HPS24\HPS24_230621_sess4';
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS25_050721_sess4 MK801 
bpath = 'F:\data\HPS25\HPS25_050721_sess4'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',17,'thetaChannel',27);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% HPS25_240621_sess2 KETAMINE
bpath = 'F:\data\HPS25\HPS25_240621_sess2'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');
    
% IPO 430
createProbe('excel_file','electrodes_coordinates_Qtrode-32ch-IPO430');
% IPO430_111021_sess1 VEHICLE
bpath = 'F:\data\IPO430\IPO430_111021_sess1'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',15,'SWChannel',[],'thetaChannel',15);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO430_131021_sess2
bpath = 'F:\data\IPO430\IPO430_131021_sess2'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO430_151021_sess3
bpath = 'F:\data\IPO430\IPO430_151021_sess3'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',17,'SWChannel',13,'thetaChannel',17);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IOP429_111021_sess1
bpath = 'F:\data\IPO429\IPO429_111021_sess1'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',13,'SWChannel',[],'thetaChannel',13);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO429_131021_sess2
bpath = 'F:\data\IPO429\IPO429_131021_sess2'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',15,'SWChannel',[],'thetaChannel',15);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO429_151021_sess3
bpath = 'F:\data\IPO429\IPO429_151021_sess3'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',15,'SWChannel',13,'thetaChannel',15);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO149
createProbe('excel_file','electrodes_coordinates_Tetrodes-16ch(4t-4c)-IPO149');
% IPO149_231021_sess3
bpath = 'F:\data\IPO149\IPO149_231021_sess3'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',6,'SWChannel',[],'thetaChannel',6);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO149_251021_sess4 KETAMINE
bpath = 'F:\data\IPO149\IPO149_251021_sess4'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',6,'SWChannel',[],'thetaChannel',6);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO149_281021_sess5
bpath = 'F:\data\IPO149\IPO149_281021_sess5'; 
preprocessSession_pablo('basepath','F:\data\IPO149\IPO149_281021_sess5','analysisPath',[],'cleanArtifacts',({[],[]}),'analogChannelsList',[],'digitalChannelsList',[],'anymaze_ttl_channel',2)
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',6,'SWChannel',[],'thetaChannel',6);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO150_091121_sess4 VEHICLE
bpath = 'F:\data\IPO150\IPO150_091121_sess4'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',2,'SWChannel',[],'thetaChannel',2);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO150_261021_sess2 MK801
bpath = 'F:\data\IPO150\IPO150_261021_sess2'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',2,'SWChannel',[],'thetaChannel',2);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO150_291021_sess3 KETAMINE
bpath = 'F:\data\IPO150\IPO150_291021_sess3'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',2,'SWChannel',2,'thetaChannel',2);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO447_181121_sess1 VEHICLE
bpath = 'F:\data\IPO447\IPO447_181121_sess1'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',12,'SWChannel',[],'thetaChannel',12);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO447_181121_sess4 KETAMINE
bpath = 'F:\data\IPO447\IPO447_241121_sess4'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',12,'SWChannel',[],'thetaChannel',12);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO447_191121_sess2 MK801
bpath = 'F:\data\IPO447\IPO447_191121_sess2'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',12,'SWChannel',12,'thetaChannel',12);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% IPO136_280122_sess1 MK801
bpath = 'F:\data\IPO136\IPO136_280122_sess1'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',30,'SWChannel',[],'thetaChannel',30);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% fCr1_220402_sess26
bpath = 'F:\data\fCr1\fCr1_220402_sess26'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'rippleChannel',56,'SWChannel',43);
processSession('basepath',bpath,'promt_hippo_layers',false);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% fCr1_220404_sess27
bpath = 'F:\data\fCr1\fCr1_220404_sess27'; 
processSession_pablo('basepath',bpath,'project','MK801Project','promt_hippo_layers',true,'rippleChannel',56,'SWChannel',43);
processSession('basepath',bpath,'promt_hippo_layers',false);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

% fCamk7_220509_sess27 VEHICLE
bpath = 'F:\data\fCamk7\fCamk7_220509_sess27'; 
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'rippleChannel',70,'thetaChannel',70);
indexNewSession_pablo('basepath',bpath,'project','MK801Project');

%% IPO197
basepath = 'F:\data\IPO197';
createDACQFiles('basepath',basepath);
%% HM PROJECT
createNSFiles('basepath','D:\HM');
bpath = 'F:\data\HM\HM_270722_sess2'; 
processSession_HM('basepath',bpath,'project','HMProject');
indexNewSession_pablo('basepath',bpath,'project','HMProject');

