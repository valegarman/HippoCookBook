
%% BatchScript_analysis_example
% place your code to run an analysis across all sessions

targetProject= 'All';

clear; close all
HCB_directory = what('HippoCookBook'); 
load([HCB_directory.path filesep 'indexedSessions.mat']);

targetProject= 'InterneuronsLibrary';
sessionNames = fieldnames(allSessions);
for ii = 1:length(sessionNames)
    sessions.basepaths{ii} = [database_path filesep allSessions.(sessionNames{ii}).path];
    sessions.project{ii} = allSessions.(sessionNames{ii}).project;
end

% EXPLORE SESSIONS AND CHECK CELL CLASSIFICATION (SAVE CHANGES AT THE END)
if ~strcmpi(targetProject,'all')
    sessions.basepaths = sessions.basepaths(strcmpi(sessions.project, targetProject));
    sessions.project = sessions.project(strcmpi(sessions.project, targetProject));
end

for ii = 1:length(sessions.basepaths)
    fprintf(' > %3.i/%3.i session \n',ii, length(sessions.basepaths)); %\n
    cd(sessions.basepaths{ii});
    try
        cd(sessions.basepaths{ii});
        
        %%% your code goes here...
        spatialModulation = getSpatialModulation('force',true);
        %%%
%         speedCorr = getSpeedCorr(pwd,'numQuantiles',20);
        close all;
    catch
        warning('Analysis was not possible!');
    end
end