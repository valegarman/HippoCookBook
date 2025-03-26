%% BatchScript_analysis_ripplesPerSubsession

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'Bibliocampus';

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);

        % ACGPeak SubSessions
        clear basepath session basename cell_metrics spikes
        basepath = pwd;
        session = loadSession;
        basename = basenameFromBasepath(basepath);
        cell_metrics = loadCellMetrics;
        parameters.forceReload = true;
        spikes = loadSpikes;
        
        if strcmpi(session.extracellular.chanCoords.layout, 'A5x12-16-Buz-lin-5mm-100-200-160-177')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
        else 
            keyboard;   
        end
        deepSuperficialfromRipple = classification_DeepSuperficial(session);
        % improving detection
        % deepSuperficialfromRipple.channelClass([18 48 17 47]) = {'Superficial'};
        cell_metrics.general.SWR = deepSuperficialfromRipple;
        deepSuperficial_ChDistance = deepSuperficialfromRipple.channelDistance; %
        deepSuperficial_ChClass = deepSuperficialfromRipple.channelClass;% cell_deep_superficial
        cell_metrics.general.deepSuperficial_file = deepSuperficial_file;
        for j = 1:cell_metrics.general.cellCount
            cell_metrics.deepSuperficial(j) = deepSuperficial_ChClass(spikes.maxWaveformCh1(j)); % cell_deep_superficial OK
            cell_metrics.deepSuperficialDistance(j) = deepSuperficial_ChDistance(spikes.maxWaveformCh1(j)); % cell_deep_superficial_distance
        end
        % save cell_metrics
        save([basename '.cell_metrics.cellinfo.mat'], "cell_metrics");
            
    end
    close all;
end