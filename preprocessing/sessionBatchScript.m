%% sessionBatchScript

%1% Transfer files and organize session's folder
updateExpFolder({'V:\data\fCck1'},'E:\data\fCck1');

%2%  Preprocessing
preprocessSession('basepath','E:\data\fCck1\fCck1_220201_sess2','getPos',false,'analogCh',1);
