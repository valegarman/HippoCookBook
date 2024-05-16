%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'All';

HCB_directory = what('HippoCookBook'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{ii}]));
        try
        
            %%% your code goes here...
            % targetFile = dir('*thetaRun_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
            % thetaChannel = thetaRunMod.detectorParams.channels;
            % targetFile = dir('*ripple_120-200.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
            % rippleChannel = rippleMod.detectorParams.channels;
            % targetFile = dir('*SW_2-10.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
            % if ~isempty(SWMod)
            %     SWChannel = SWMod.detectorParams.channels;
            % else
            %     SWChannel = [];
            % end
            % 
            % computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'hgammaChannel',thetaChannel,'lgammaChannel',thetaChannel);
            % clear rippleChannel SWChannel thetaChannel
            getAverageCCG('force',true);
            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end

%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions