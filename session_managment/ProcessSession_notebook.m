
%% ProcessSession_notebook
% Andrea Gallardo and Manu Valero, 2023

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
ignoretime = [];
ThetaChannels = [20 46 19 45 21 43]; % choose theta channels (ideally SLM)
SWChannels = [30 36 26]; % choose slow wave channels (ideally superficial cortex)

% if EMG is not been quantified correctly, try discarting channels with
bz_EMGFromLFP(pwd,'rejectChannels',[51 58 52 57 49 60 50 59 54 55 53 56],'overwrite', true);

targetFile = (dir('*optogeneticPulses.events.mat'));
if ~isempty(targetFile)
    pulses = importdata(targetFile.name);
else
    pulses.stimulationEpochs = [];
end
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[pulses.stimulationEpochs; ignoretime], 'overwrite', true, 'ThetaChannels', ThetaChannels, 'SWChannels', SWChannels,'rejectChannels',[51 58 52 57 49 60 50 59 54 55 53 56]);

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
rippleChannel = [4];
SWChannel = [44];
eventSpikeThreshold_shanks = [4 5]; % which shanks will be accounted for the spike threshold 
rippleMasterDetector_threshold = [1.5 3.5]; % [1.5 3.5]
eventSpikeThreshold = 1;
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'removeOptogeneticStimulation',false,'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold_shanks', eventSpikeThreshold_shanks,'eventSpikeThreshold',eventSpikeThreshold);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
getSpikesRank('events','ripples');

%% 5. Theta epochs
% Revise channel definition, theta band in thetaEpochs.png and cells
% rhytmicity. If bad, you can change useCSD_for_theta_detection to false,
% or change powerThreshold, even the channel
useCSD_for_theta_detection = true;
powerThreshold = 1;
channel = [18];
thetaEpochs = detectThetaEpochs('force',true,'useCSD',useCSD_for_theta_detection,'powerThreshold',powerThreshold,'channel', channel);

%% 6. Phase modulation
% NOTE!! If you have re-detected ripples or theta, you must run this code
% again with the same channel definition!!!
thetaChannel = [15];
[phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel);
computeCofiringModulation;

%% 7. Brain region
% Revise assignBrainRegion output. If disagreements or out of date, run again
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples'); % hfo slowOscilations [-.5 .5]

%% 8. Cell metrics
% Basically, if you have change anything above, a good practice is running
% this code again (but unless you merge or discard cells in phy, you don´t
% have to validate the monosynpatic connections again)
    
if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
    load(file.name);
end
excludeManipulationIntervals = optoPulses.stimulationEpochs;

cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true);

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
getSessionArmChoice('task','alternation','leftArmTtl_channel',2,'rightArmTtl_channel',3,'homeDelayTtl_channel',4);
behaviour = getSessionLinearize('forceReload',true);  
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

% -------------------------------------------------------------------------



