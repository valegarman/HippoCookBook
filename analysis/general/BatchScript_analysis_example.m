
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
        behaviour = getSessionLinearize;
        psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',15,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01);
        psth_rReward = spikesPsth([behaviour.events.rReward],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',15,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01);
        psth_intersection = spikesPsth([behaviour.events.intersection],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',15,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01);
        psth_startPoint = spikesPsth([behaviour.events.startPoint],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',15,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01);

        behaviour.psth_lReward = psth_lReward;
        behaviour.psth_rReward = psth_rReward;
        behaviour.psth_intersection = psth_intersection;
        behaviour.psth_startPoint = psth_startPoint; 
        behavior = behaviour; % british to american :)
        save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior');
        clear behavior behaviour
    catch
        warning('Analysis was not possible!');
    end
end