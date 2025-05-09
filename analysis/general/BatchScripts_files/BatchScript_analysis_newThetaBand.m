%% BatchScript_analysis_ripplesPerSubsession

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'Bibliocampus';

for ii = 1:length(sessionsTable.SessionName)
    if contains(sessionsTable.Project(ii), targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}]);

        % 
        clear basepath
        [phaseMod] = computePhaseModulation('saveMat', false);

        thetaMod = phaseMod.theta;
        save([basenameFromBasepath(pwd) '.theta_5-12.PhaseLockingData.cellinfo.mat'], 'thetaMod');

        thetaREMMod = phaseMod.thetaREMMod;
        save([basenameFromBasepath(pwd) '.thetaREM_5-12.PhaseLockingData.cellinfo.mat'], 'thetaREMMod');

        thetaRunMod = phaseMod.thetaRunMod;
        save([basenameFromBasepath(pwd) '.thetaRun_5-12.PhaseLockingData.cellinfo.mat'], 'thetaRunMod');
    end
    close all;
end