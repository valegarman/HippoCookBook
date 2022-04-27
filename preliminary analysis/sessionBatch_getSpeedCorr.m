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
    
%     speedCorrs = getSpeedCorr(basepath);
end
    
    
    stop_indx = find(speed.velocity < 0.5);
    start_indx = find(speed.velocity > 5);
    stop_start = zeros(length(speed.velocity),1);
    stop_start(stop_indx) = -1;
    stop_start(start_indx) = 1;
    
    
    stop_indx = find(speed.acceleration < 0.5);
    start_indx = find(speed.acceleration > 5);
    stop_start = zeros(length(speed.acceleration),1);
    stop_start(stop_indx) = -1;
    stop_start(start_indx) = 1;
    d = [0; diff(stop_start)];
    timestamps_1 = find(d == 1);
    timestamps_2 = find(d == -1);
    timestamps_3 = find(d == 2);
    
    spikesPsth(timestamps_1,'winSize',);
    
    
end