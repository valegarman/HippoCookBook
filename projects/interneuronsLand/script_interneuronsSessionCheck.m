%% Steps to improve session

%1. Assing probe
plotProbe('force',true);

%2. check hippocampal layers, if not good run:
[hippocampalLayers] = getHippocampalLayers('force',true,'promt',true);

%3. Loot at ripples, if not good run:
SWChannel = []; 
rippleChannel = [];
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'useCSD',false);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);

%4. assing layers
session = assignBrainRegion;

%5. get ACG peak
getACGPeak('force',true);

%6. Loot at thetaEpochs, if not good run:
SWChannel = []; 
rippleChannel = [];
[phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel);
computeCofiringModulation;

%7. Speed correla
speedCorr = getSpeedCorr('force',true);

%8. Re-ran processCellMetrics
session = loadSession;
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end

cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,...
    'forceReloadWaveformBasedMetrics',false);

%9. Make summary
plotSummary;

%% Steps if cell is cleanned on phy
spikes = loadSpikes('forceReload',true);
spikeFeatures;
optogeneticResponses = getOptogeneticResponse('numRep',500,'force',true);

psthUD = spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true);

psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);

[phaseMod] = computePhaseModulation;
computeCofiringModulation;

behaviour = getSessionLinearize;
spikes = loadSpikes;
firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
firingTrialsMap = firingMapPerTrial('force',true);
spatialModulation = getSpatialModulation('force',true);

behaviour = getSessionLinearize;
psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',100,'saveMat',false,...
    'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
psth_rReward = spikesPsth([behaviour.events.rReward],'numRep',100,'saveMat',false,...
    'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
psth_reward = spikesPsth([behaviour.events.lReward; behaviour.events.rReward],'numRep',100,'saveMat',false,...
    'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);

psth_intersection = spikesPsth([behaviour.events.intersection],'numRep',100,'saveMat',false,...
    'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
psth_startPoint = spikesPsth([behaviour.events.startPoint],'numRep',100,'saveMat',false,...
    'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);

behaviour.psth_lReward = psth_lReward;
behaviour.psth_rReward = psth_rReward;
behaviour.psth_reward = psth_reward;
behaviour.psth_intersection = psth_intersection;
behaviour.psth_startPoint = psth_startPoint; 
behavior = behaviour; % british to american :)
save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior');
    
speedCorr = getSpeedCorr('numQuantiles',20,'force',true);

getSpikesReturnPlot('force',true);

session = loadSession;
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end

cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,...
    'forceReloadWaveformBasedMetrics',false);

getAverageCCG('force',true);

getACGPeak('force',true);

plotSummary;


%% optionals
session = sessionTemplate(pwd,'showGUI',true);


spikeFeatures;

spikes = loadSpikes('forceReload',true);


%
excludeManipulationIntervals = [0 Inf];
excludeManipulationIntervals = SubtractIntervals(excludeManipulationIntervals,[1.2E4 2E4])
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',false,...
    'forceReloadWaveformBasedMetrics',true);

spikeFeatures('skipStimulationPeriods', true);

file = dir([basenameFromBasepath(pwd),'.optogeneticPulses.events.mat']); load(file.name);
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',optoPulses.stimulationEpochs, 'overwrite', true);
TheStateEditor(session.general.name);
bz_ThetaStates(pwd);