%% ======== CHECKING INTERNEURONS SESSIONS ================
% GENERAL ANALYSIS
% 1. He incluido en getSummaryPerCell phase modulation de theta durante Run
% y REM. Creo que visualmente puede quedar mejor sobre todo cuando vemos
% bimodalidades ( por ejemplo en fNkx9_sess9)
%% 'fnkx9_200818_sess2'

basepath = 'Z:\data\fNkx9\fNkx9_200818_sess2';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getSummaryPerCell;


%%
basepath = 'Z:\data\fNkx8\fNkx8_200817_sess1';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getSummaryPerCell;

%%
basepath = 'Z:\data\fNkx8\fNkx8_200819_sess3';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getSummaryPerCell;

%% 'fNkx8_200826_sess8'
basepath = 'Z:\data\fNkx8\fNkx8_200826_sess8';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();

%% 'fNkx9_200825_sess7' Responsive cells that look like pyr ( even low firing rate)
basepath = 'Z:\data\fNkx9\fNkx9_200825_sess7';
cd(basepath);
% Brain Regions
session = assignBrainRegion();
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getACGPeak;
getSummaryPerCell;


%% 'fNkx9_200827_sess9' Only 1 responsive cell. Check speed relationship, Check theta phase (bimodal), REM vs RUN
basepath = 'Z:\data\fNkx9\fNkx9_200827_sess9';
cd(basepath);
% Channel 17 is not a bad Channel
session = sessionTemplate(basepath,'showGUI',true);
indexNewSession('basepath',basepath,'force_loadingSpikes', false,'force_analogPulsesDetection',false,'promt_hippo_layers',true,'analogChannelsList',65,'removeDatFiles',false,'indexing',false);

% Brain Regions
session = assignBrainRegion();
% Re run getSpeedCorr
speedCorr = getSpeedCorr('numQuantiles',20,'force',true);
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getACGPeak;
getSummaryPerCell;

%% 'fNkx9_200909_sess17' 
basepath = 'Z:\data\fNkx9\fNkx9_200909_sess17';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Channel 17 is not a bad Channel
session = sessionTemplate(basepath,'showGUI',true);
indexNewSession('basepath',basepath,'force_loadingSpikes', false,'force_analogPulsesDetection',false,'promt_hippo_layers',true,'analogChannelsList',65,'removeDatFiles',false,'indexing',false);
% Brain Regions
session = assignBrainRegion();
% Re run getSpeedCorr
speedCorr = getSpeedCorr('numQuantiles',20,'force',true);
% Re run cell_metrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
% Get summary per cell
getACGPeak;
getSummaryPerCell;