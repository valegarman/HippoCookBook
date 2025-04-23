%% BatchScript_analysis_ripplesPerSubsession

clear; close all
HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions
targetProject = 'Bibliocampus';

for ii = 214:length(sessionsTable.SessionName)
    if contains(sessionsTable.Project(ii), targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}]);

        % 
        clear basepath session basename cell_metrics spikes
        basepath = pwd;
        session = loadSession;
        basename = basenameFromBasepath(basepath);
        cell_metrics = loadCellMetrics;
        parameters.forceReload = true;
        spikes = loadSpikes;
        deepSuperficial_file = fullfile(basepath, [basename,'.deepSuperficialfromRipple.channelinfo.mat']);
        
        if strcmpi(session.extracellular.chanCoords.layout, 'A5x12-16-Buz-lin-5mm-100-200-160-177')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'A5x12-16-Buz-lin-5mm-100-200-160-177';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'CambridgeNeurotech-E1-64ch(4shank-ege)')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'CambridgeNeurotech-E1-64ch(4shank-ege)';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'uLED-12LED-32Ch-4Shanks')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'uLED-12LED-32Ch-4Shanks';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'DiagnosticBiochip-128-6-128ch')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'DiagnosticBiochip-128-6-128ch';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'CambridgeNeurotech-H2-64ch(2shank-ege)')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'CambridgeNeurotech-H2-64ch(2shank-ege)';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'poly2')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'A5x12-16-Buz-lin-5mm-100-200-160-177';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'staggered')
            session.extracellular.chanCoords.probe = session.animal.probeImplants{1}.probe;
        elseif strcmpi(session.extracellular.chanCoords.layout, 'CambridgeNeurotech-H3-64ch')
            session.extracellular.chanCoords.layout = 'linear';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'CambridgeNeurotech-H3-64ch';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'A3x8-16-Buz-lin-5mm-50-150-160-703')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'A3x8-16-Buz-lin-5mm-50-150-160-703';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'DiagnosticBiochip-128-6-128ch&uLED_12LED-32Ch-4Shanks')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'DiagnosticBiochip-128-6-128ch&uLED_12LED-32Ch-4Shanks';
        elseif strcmpi(session.extracellular.chanCoords.layout, 'A1x32-Poly3-5mm-25s-177')
            session.extracellular.chanCoords.layout = 'staggered';
            session.extracellular.chanCoords.verticalSpacing = 20;
            session.extracellular.chanCoords.probe = 'A1x32-Poly3-5mm-25s-177';
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
        save([basename '.session.mat'], "session");
       
        computeRippleReversal;
    end
    close all;
end