%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project
clear; close all

list_of_paths = {'Y:\fCamk1\fCamk1_200827_sess9', 'Y:\fCamk1\fCamk1_200901_sess12', 'Y:\fCamk1\fCamk1_200902_sess13',...
    'Y:\fCamk1\fCamk1_200904_sess15','Y:\fCamk1\fCamk1_200908_sess16','Y:\fCamk1\fCamk1_200909_sess17','Y:\fCamk1\fCamk1_200910_sess18',...
    'Y:\fCamk1\fCamk1_200911_sess19','Y:\fCamk3\fCamk3_201028_sess10_cleanned','Y:\fCamk3\fCamk3_201029_sess11_cleanned',...
    'Y:\fCamk3\fCamk3_201030_sess12','Y:\fCamk3\fCamk3_201102_sess13','Y:\fCamk3\fCamk3_201103_sess14','Y:\fCamk3\fCamk3_201111_sess20',...
    'Y:\fCamk3\fCamk3_201113_sess22','Y:\fCamk3\fCamk3_201105_sess16','Y:\fCamk3\fCamk3_201106_sess17','Y:\fCamk3\fCamk3_201110_sess19',...
    'Y:\fCamk3\fCamk3_201109_sess18', 'Y:\fCamk5\fCamk5_210406_sess10', 'Y:\fCamk5\fCamk5_210408_sess12', ...
    'Y:\fCamk5\fCamk5_210412_sess14', 'Y:\fCamk5\fCamk5_210415_sess17'};

for ii = 1:length(list_of_paths)
    fprintf(' > %3.i/%3.i session \n', ii, length(list_of_paths));
    cd(adapt_filesep(list_of_paths{ii}));
    try
        %%% your code goes here...
        spatialModulation = getSpatialModulation('force',true);
        %%%
        close all;
    catch
        warning('Analysis was not possible for session: %s', list_of_paths{ii});
    end
end
