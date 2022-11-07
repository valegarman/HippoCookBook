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


%% SPATIAL SUBICULUM PROJECT
createProbe('excel_file','electrodes_coordinates_Buzsaki64(64 ch, 8 shanks, staggered)');
bpath = 'F:\data\HPS22\HPS22_210521_sess17';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',39,'SWChannel',13,'thetaChannel',39,'excludeAnalysis',{'getHippocampalLayers'});
indexNewSession_pablo('basepath',bpath,'project','Subiculum Project');


%% MK801 Project
% HPS22_100621_sess26
bpath = 'F:\data\HPS22\HPS22_100621_sess26';
% processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',24,'SWChannel',13,'thetaChannel',24,'excludeAnalysis',{'getHippocampalLayers'});
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',24,'SWChannel',13,'thetaChannel',24);
indexNewSession_pablo('basepath',bpath,'project','MK801 Project');





