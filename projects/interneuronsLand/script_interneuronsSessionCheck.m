%% Steps to improve session

%1. Assing probe
plotProbe('force',true);

%2. check hippocampal layers, if not good run:
[hippocampalLayers] = getHippocampalLayers('force',true,'promt',true);

%3. Loot at ripples, if not good run:
SWChannel = []; 
rippleChannel = [];
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'useCSD',false,'threshold',[1.25 3]);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);

%4. assing layers
session = assignBrainRegion;

%5. get ACG peak
getACGPeak('force',true);

%6. Loot at thetaEpochs, if not good run:
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false);

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
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',optoPulses.stimulationEpochs, 'overwrite', true,'ThetaChannels',[20 46 56],'SWChannels',[30 36 26]);
TheStateEditor(session.general.name);
bz_ThetaStates(pwd);


% to load uLED sessions

session = loadSession(basepath);
if ~isfield(session, 'analysisTags')
    session.analysisTags = [];
end
session.analysisTags.digital_optogenetic_channels = [11 12 13 14 15 16];
session.analysisTags.analog_optogenetic_channels = [3 4 5 6 7 8];
session.analysisTags.bazler_ttl_channel = 1;
save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');

session = sessionTemplate(pwd,'showGUI',true); % analog
plotProbe('force',true); % choose probe

spikes = loadSpikes('forceReload',true);

optogeneticResponses = getOptogeneticResponse('numRep',500,'force',true); % this takes literally forever...
spikeFeatures;
bz_ThetaStates(pwd);

powerProfile_theta = powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true);
powerProfile_gamma = powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true);
powerProfile_hfo = powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true);

[hippocampalLayers] = getHippocampalLayers('force',true,'promt',true);

ripples = rippleMasterDetector('rippleChannel',[],'SWChannel',[],'force',true,'removeRipplesStimulation',false);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false);
[phaseMod] = computePhaseModulation('rippleChannel',[],'SWChannel',[],'skipStimulationPeriods',true);
computeCofiringModulation;
session = assignBrainRegion;

try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
            file = dir([session.general.name,'.optogeneticPulses.events.mat']);
            load(file.name);
    end
            excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
        warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true);

getACGPeak('force',true,'lightPulseDuration',0.02);

getAverageCCG('force',true);

getSpikesReturnPlot('force',true);

firingTrialsMap = firingMapPerTrial('force',true);
spatialModulation = getSpatialModulation('force',true);
behaviour = getSessionLinearize;
psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
psth_rReward = spikesPsth([behaviour.events.rReward],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
psth_reward = spikesPsth([behaviour.events.lReward; behaviour.events.rReward],'numRep',100,'saveMat',false,...
            'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
        
behaviour.psth_lReward = psth_lReward;
behaviour.psth_rReward = psth_rReward;
behaviour.psth_reward = psth_reward;
behaviour.psth_intersection = NaN;
behaviour.psth_startPoint = NaN; 
behavior = behaviour; % british to american :)
save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior');

speedCorr = getSpeedCorr('numQuantiles',20,'force',true);

plotSummary('cells_responsive_to_any_pulse',true);

 
   



