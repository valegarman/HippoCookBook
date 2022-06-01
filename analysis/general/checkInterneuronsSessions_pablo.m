%% ======== CHECKING INTERNEURONS SESSIONS ================
%%
basepath = 'Z:\data\fNkx9\fNkx9_200818_sess2';
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
getSummaryPerCell;


%%
basepath = 'Z:\data\fNkx8\fNkx8_200817_sess1';
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
getSummaryPerCell;

%%
basepath = 'Z:\data\fNkx8\fNkx8_200819_sess3';
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
getSummaryPerCell;

%%
basepath = 'Z:\data\fNkx8\fNkx8_200826_sess8';
cd(basepath);
% Brain Regions
session = assignBrainRegion();



