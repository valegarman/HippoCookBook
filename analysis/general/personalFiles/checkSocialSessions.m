%% ======== CHECKING SOCIAL SESSIONS ================

% Change brain regions to CA2 or CA3 and rerun cell_metrics

%% FL4_080322_sess1

bpath = 'D:\FLR\FL4\FL4_080322_sess1';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
cell_metrics = CellExplorer('basepath',pwd);
getAverageCCG('force',true);
close all;

%% FL4_090322_sess2

bpath = 'D:\FLR\FL4\FL4_090322_sess2';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
cell_metrics = CellExplorer('basepath',pwd);
getAverageCCG('force',true);
close all;

%% FL4_100322_sess3

bpath = 'D:\FLR\FL4\FL4_100322_sess3';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
cell_metrics = CellExplorer('basepath',pwd);
getAverageCCG('force',true);
close all;

%% FL4_110322_sess4

bpath = 'D:\FLR\FL4\FL4_110322_sess4';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
cell_metrics = CellExplorer('basepath',pwd);
getAverageCCG('force',true);
close all;

%% FL4_170322_sess5

bpath = 'D:\FLR\FL4\FL4_170322_sess5';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);
close all;

%% FL4_180322_sess6

bpath = 'D:\FLR\FL4\FL4_180322_sess6';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);
close all;

%% FL4_210322_sess7

bpath = 'D:\FLR\FL4\FL4_210322_sess7';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);
close all;

%% FL4_220322_sess8

bpath = 'D:\FLR\FL4\FL4_220322_sess8';
cd(bpath);

session = loadSession();
% session.analysisTags.thetaChannel = 27;
% save([session.general.name,'.session.mat'],'session');
% 
% % 8.3 Theta intervals
% thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',27);
% % LFP-spikes modulation
% [phaseMod] = computePhaseModulation('rippleChannel',27,'SWChannel',32,'thetaChannel',27,'lgammaChannel',27,'hgammaChannel',27);    

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);
close all;

%% FL3_020322_sess1

bpath = 'D:\FLR\FL3\FL3_020322_sess1';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);


%% FL3_030322_sess2

bpath = 'D:\FLR\FL3\FL3_030322_sess2';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);


%% FL3_100322_sess3

bpath = 'D:\FLR\FL3\FL3_100322_sess3';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);

%% FL3_110322_sess4

bpath = 'D:\FLR\FL3\FL3_110322_sess4';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);


%% FL3_160322_sess5

bpath = 'D:\FLR\FL3\FL3_160322_sess5';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);

%% FL3_170322_sess6

bpath = 'D:\FLR\FL3\FL3_170322_sess6';
cd(bpath);

session = loadSession();

session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    excludeManipulationIntervals = [];
end
cell_metrics = ProcessCellMetrics('session', session,'manualAdjustMonoSyn',false,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);
getAverageCCG('force',true);

