%% BatchScript_analysis

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'Bibliocampus';

for ii = 220:length(sessionsTable.SessionName)
    %if contains(sessionsTable.Project(ii), targetProject)% || strcmpi('all', targetProject)
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

        lgamma = phaseMod.lgamma;
        save([basenameFromBasepath(pwd) '.lgamma_20-50.PhaseLockingData.cellinfo.mat'], 'lgamma');

        hgamma = phaseMod.hgamma;
        save([basenameFromBasepath(pwd) '.hgamma_50-100.PhaseLockingData.cellinfo.mat'], 'hgamma');

        fgamma = phaseMod.fgamma;
        save([basenameFromBasepath(pwd) '.fgamma_100-140.PhaseLockingData.cellinfo.mat'], 'fgamma');
    %end
    close all;
end