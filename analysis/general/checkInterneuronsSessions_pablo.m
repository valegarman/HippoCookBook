%% ======== CHECKING INTERNEURONS SESSIONS ================
% GENERAL ANALYSIS
% 1. He incluido en getSummaryPerCell phase modulation de theta durante Run
% y REM. Creo que visualmente puede quedar mejor sobre todo cuando vemos
% bimodalidades ( por ejemplo en fNkx9_sess9)


%% fNkx9_200820_sess4 FALTABA POR INCLUIR ESTA SESION EN EL CSV
basepath = 'Z:\data\fNkx9\fNkx9_200820_sess4';
cd(basepath);
bpath = 'Z:\data\fNkx9\fNkx9_200820_sess4';
indexNewSession('basepath',bpath,'force_loadingSpikes',true,'analogChannelsList',65);
% 8.3 Theta intervals
thetaEpochs = detectThetaEpochs('force',true);
% LFP-spikes modulation
rippleChannel = [];
SWChannel = [];
[phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel);
computeCofiringModulation;
% 8.2 Ripples
rippleChannel = [];
SWChannel = [];
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);
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
% Summary
getAverageCCG('force',true);
getACGPeak;
getSummaryPerCell;

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
%% fNkx8_200817_sess1
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

%% fNkx8_200819_sess3
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
% getHippocampal Layers
promt_hippo_layers = true;
[hippocampalLayers] = getHippocampalLayers('force',true,'promt',promt_hippo_layers);
% 8.2 Ripples
ripples = rippleMasterDetector('rippleChannel',[],'SWChannel',44,'force',true);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);
% 8.3 Theta intervals
thetaEpochs = detectThetaEpochs('force',true);
% 9. Phase Modulation
% LFP-spikes modulation
[phaseMod] = computePhaseModulation('rippleChannel',[],'SWChannel',44);
computeCofiringModulation;
% Brain Regions
session = assignBrainRegion();
% 10. Cell metrics
% Exclude manipulation intervals for computing CellMetrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'excludeMetrics',{'deepSuperficial'});
speedCorr = getSpeedCorr('numQuantiles',20,'force',true);
getACGPeak('force',true);
getSummaryPerCell;
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
session = loadSession(basepath);
% Channel 17 is not a bad Channel
session = sessionTemplate(basepath,'showGUI',true);
indexNewSession('basepath',basepath,'force_loadingSpikes', false,'force_analogPulsesDetection',false,'promt_hippo_layers',true,'analogChannelsList',65,'removeDatFiles',false,'indexing',false);

% Brain Regions
session = assignBrainRegion();
% Re run getSpeedCorr
speedCorr = getSpeedCorr('numQuantiles',20,'force',true);
% 8.3 Theta intervals
thetaEpochs = detectThetaEpochs('force',true);
% 9. Phase Modulation
% LFP-spikes modulation
[phaseMod] = computePhaseModulation('rippleChannel',[],'SWChannel',[]);
computeCofiringModulation;
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

%% fNkx8_200921_sess12. Problems with phaseModulation with pyr ( is a line)
basepath = 'Z:\data\fNkx8\fNkx8_200901_sess12';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();
% LFP-spikes modulation
[rippleMod,SWMod,thetaMod,lgammaMod,hgammaMod,thetaRunMod,thetaREMMod] = computePhaseModulation('rippleChannel',[],'SWChannel',[]);
computeCofiringModulation;
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

%% fNkx8_200902_sess13
basepath = 'Z:\data\fNkx8\fNkx8_200902_sess13';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
% Brain Regions
session = assignBrainRegion();
bpath = 'Z:\data\fNkx8\fNkx8_200902_sess13';
indexNewSession('basepath',bpath,'promt_hippo_layers',true,'analogChannelsList',65,'bazler_ttl_channel',1);

% Get summary per cell
getACGPeak;
getSummaryPerCell;

bpath = 'Z:\data\fNkx8\fNkx8_200902_sess13';
indexNewSession('basepath',bpath,'promt_hippo_layers',true,'analogChannelsList',65);
%% fNkx9_200821_sess5
basepath = 'Z:\data\fNkx9\fNkx9_200821_sess5';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx11_201029_sess11
basepath = 'Z:\data\fNkx11\fNkx11_201029_sess11';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx8_200908_sess16
basepath = 'Z:\data\fNkx8\fNkx8_200908_sess16';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx11_201104_sess15
basepath = 'Z:\data\fNkx11\fNkx11_201104_sess15';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx11_201015_sess1
basepath = 'Z:\data\fNkx11\fNkx11_201015_sess1';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fSst1_190603_sess33
basepath = 'Z:\data\fSst1\fSst1_190603_sess33';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv3_210305_sess5
basepath = 'Z:\data\fPv3\fPV3_210305_sess5';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv3_210304_sess4
basepath = 'Z:\data\fPv3\fPV3_210304_sess4';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx8_200915_sess21
basepath = 'Z:\data\fNkx8\fNkx8_200915_sess21';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv4_210310_sess3
basepath = 'Z:\data\fPv4\fPv4_210310_sess3';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv4_210309_sess2
basepath = 'Z:\data\fPv4\fPv4_210309_sess2';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fNkx11_201019_sess3
basepath = 'Z:\data\fNkx11\fNkx11_201019_sess3';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv3_210302_sess2
basepath = 'Z:\data\fPv3\fPV3_210302_sess2';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);

%% fPv4_210311_sess4 
basepath = 'Z:\data\fPv4\fPv4_210311_sess4';
cd(basepath);
% The State Editor
session = loadSession(basepath);
TheStateEditor_temp(session.general.name);
