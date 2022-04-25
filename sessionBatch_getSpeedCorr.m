% loadProjectResults('project', 'InterneuronsLibrary');
    
clear; close all
HCB_directory = what('HippoCookBook');
load([HCB_directory.path filesep 'indexedSessions.mat']);

targetProject= 'InterneuronsLibrary';
sessionNames = fieldnames(allSessions);

for i = 1:length(sessionNames)
    cd(['Z:\data', filesep, allSessions.(sessionNames{i}).path]);
    basepath = pwd;
    session = sessionTemplate(basepath);
    if ~isempty(dir([session.general.name, '.Speed.events.mat']))
        forceReload = true;
    else
        forceReload = false;
    end
    getSessionTracking('convFact',0.1149,'roiTracking','manual','forceReload',forceReload);
    
    speedCorrs = getSpeedCorr(basepath);
    
    getSpeedAccelerationOnset()
    
    
    
    
end