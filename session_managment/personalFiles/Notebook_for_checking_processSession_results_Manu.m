
%% Notebook_for_checking_processSession_results
% Andrea Gallardo and Manu Valero, 2023
%% 0. Check metadata
gui_session;

%% 1. Optogenetic responses and cell identity
% Check light responses. If two neurons are tagged, but they show a similar
% waveform and same max channel, they are probably the same cell.
% The only way of fixing this is merging them in Phy!!!!!!! (and then run
% processSession again)

%% 2. Brain state scoring
% Inspect content in StateScoreFigures. Three nice clusters must be visible
% in SSCluster3D.jpg. Otherwise, try following:

% 2. 1 Edit scoring with 
session = loadSession;
TheStateEditor(session.general.name);
% press 'a', click Sticky and move thresholds for improving detections of
% bimodalities. 
% 
% 2.2 If no clear bimodalities are visible, run scoring with different
% channels and include noisy ignoretime epochs
ThetaChannels = 60; % choose theta channels (ideally SLM)
SWChannels = 36; % choose slow wave channels (ideally superficial cortex)
ignoretime = [9000 15000];
% if EMG is not been quantified correctly, try discarting channels with
% bz_EMGFromLFP(pwd,'rejectChannels',[],'overwrite', true,'ignoretime',ignoretime);

targetFile = (dir('*optogeneticPulses.events.mat'));
if ~isempty(targetFile)
    pulses = importdata(targetFile.name);
else
    pulses.stimulationEpochs = [];
end
ignoretime = [pulses.stimulationEpochs; ignoretime]; % [pulses.stimulationEpochs; ignoretime]
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',ignoretime, 'overwrite', true, 'ThetaChannels', ThetaChannels, 'SWChannels', SWChannels,'rejectChannels',[]);
% SleepScoreMaster(pwd,'noPrompts',true,'overwrite', true, 'ThetaChannels', ThetaChannels, 'SWChannels', SWChannels,'rejectChannels',[],'scoretime',[1 440*60]);
 
% 2.3 As the last resource, you can edit the epochs in TheStateEditor. IT
% IS NOT WORKING!!!!!
TheStateEditor(session.general.name);

%% 3. Hippocampal layers
% Revise hippocapmal layer definition. Upon disagrements, run again
hippocampalLayers = getHippocampalLayers('force',true,'promt',true);

%% 4. Ripples
% Revise number of ripples, shape and channel. You can also specifiy a
% different trhresold or restrict the shanks for the spikeThreshold
% analysis.

excludeIntervals = [];
rippleChannel = [];
SWChannel = [];
noiseChannel = [];
eventSpikeThreshold_shanks = [3]; % which shanks will be accounted for the spike threshold 
eventSpikeThreshold = .8; % .5
rippleMasterDetector_threshold = [1.5 3.5]; % [1.5 3.5]
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'skipStimulationPeriods',false,'thresholds',...
    rippleMasterDetector_threshold,'eventSpikeThreshold_shanks', eventSpikeThreshold_shanks,'eventSpikeThreshold',eventSpikeThreshold,'excludeIntervals',excludeIntervals,'noise',noiseChannel); 
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
% psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);
getSpikesRank('events','ripples');

%% 5. Theta epochs
% Revise channel definition, theta band in thetaEpochs.png and cells
% rhytmicity. If bad, you can change useCSD_for_theta_detection to false,
% or change powerThreshold, even the channel
channel = 45;
useCSD_for_theta_detection = true;
powerThreshold = 1.2;% .8
thetaEpochs = detectThetaEpochs('force',true,'useCSD',useCSD_for_theta_detection,'powerThreshold',powerThreshold,'channel', channel);

%% 6. Phase modulation
% NOTE!! If you have re-detected ripples or theta, you must run this code
% again with the same channel definition!!!
thetaChannel = [];
hgammaChannel = [];
lgammaChannel = [];
[phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'hgammaChannel',thetaChannel,'lgammaChannel',thetaChannel);
computeCofiringModulation;
% [phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'hgammaChannel',thetaChannel,'lgammaChannel',thetaChannel,'restrict_to_manipulation',true);

%% 7. Brain region
% Revise assignBrainRegion output. If disagreements or out of date, run again
assignBrainRegion('showPowerProfile','theta','showEvent','ripples'); % hfo slowOscilations [-.5 .5]

%% 8. Cell metrics
% Basically, if you have change anything above, a good practice is running
% this code again (but unless you merge or discard cells in phy, you don´t
% have to validate the monosynpatic connections again)
session = loadSession;
if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
    load(file.name);
end
excludeManipulationIntervals = optoPulses.stimulationEpochs;

cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'getWaveformsFromDat', true);
% ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'restrictToIntervals',...
%     [0 1.2957e+04],'manualAdjustMonoSyn',false); % uLEDResponses.restricted_interval
% cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,...
%     'restrictToIntervals',[1.2957e+04 2.0668e+04],'manualAdjustMonoSyn',false,'saveAs','cell_metrics_post'); % uLEDResponses.restricted_interval

getACGPeak('force',true);

getAverageCCG('force',true);
    
getSpikesReturnPlot('force',true);

%% 9. Spatial modulation
% Open session folder with behaviour (with a video, normaly .avi, file). If
% no folder is present, you must run this code. You can play with the
% LED_threshold level (threshold for LED detection), and check if TTL
% channel are well defined. 
LED_threshold = 0.98;
spikes = loadSpikes;
getSessionTracking('roiTracking','manual','forceReload',false,'LED_threshold',LED_threshold);
% only if the animal run a figure-eight maze behavior
%getSessionArmChoice('task','alternation','leftArmTtl_channel',2,'rightArmTtl_channel',3,'homeDelayTtl_channel',4);
behaviour = getSessionLinearize('forceReload',true,'maze','tMaze');  
firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
firingTrialsMap = firingMapPerTrial('force',true);
spatialModulation = getSpatialModulation('force',true);

% if rReward or lReward dots are not clusterized in front of the IRSensor
% position, pray (work in progress)

%% 10. Summary per cell
% Check how many cells have been included for the further analysis
% Also, if substancial changes were done before, run again!!
% (Summary_cell_1, 2 etc). If bad cells, run again. 
plotSummary('showTagCells',true);
close all;

% -------------------------------------------------------------------------



