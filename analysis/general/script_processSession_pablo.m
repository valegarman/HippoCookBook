%% SOCIAL PROJECT
createProbe('excel_file','electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-4');

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
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',27,'SWChannel',32,'thetaChannel',28);
indexNewSession_pablo('basepath',bpath,'project','Social Project');




%% SPATIAL SUBICULUM PROJECT
createProbe('excel_file','electrodes_coordinates_Buzsaki64(64 ch, 8 shanks, staggered)');
bpath = 'F:\data\HPS22\HPS22_210521_sess17';
processSession_pablo('basepath',bpath,'promt_hippo_layers',true,'anymaze_ttl_channel',2,'rippleChannel',39,'SWChannel',13,'thetaChannel',39,'excludeAnalysis',{'getHippocampalLayers'});