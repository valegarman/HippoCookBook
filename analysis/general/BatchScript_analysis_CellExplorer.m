%% BatchScript_analysis_ripplesPropertiesPerSubSession
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'MK801Project';
cd('F:\data');
database_path = 'F:\data';
HCB_directory = what('MK801Project'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_MK801Project.csv']); % the variable is called allSessions
forceReload = false;

win_resp = [-0.025 0.025];

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        
        cd([database_path filesep sessionsTable.Path{ii}]);
        
        if ~exist(dir(['\SummaryFigures\Summary_Checked.png']))
            session = loadSession();
            
            % Checking States
            TheStateEditor(session.general.name);

            % Rerun thetaEpochs and modulation
            % 8.3 Theta intervals
            file = dir(['*thetaEpochs.states.mat']); load(file.name);
            thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);
            
            % LFP-spikes modulation
            file = dir(['*ripples.events.mat']); load(file.name);
            file = dir(['*sharpwaves.events.mat']); load(file.name);
            
            [phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionhannel,'SWChannel',SWChannel);
            computeCofiringModulation;
            
             % LFP-spikes modulation per subsession
            [phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',rippleChannel,'SWChannel',SWChannel);

            % Rerun ProcessCellmetrics
            cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

            % Checking CellExplorer
            cell_metrics = loadCellMetrics('basepath',pwd);
            cell_metrics = CellExplorer('metrics',cell_metrics);

            % Plotting again
            plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
            close(gcf);
        end
        
    end
    close all;
end