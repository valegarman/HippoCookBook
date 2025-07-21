
%% Notebook_for_checking_processSession_results
% Andrea Gallardo and Manu Valero, 2023
%% 0. Check metadata
gui_session;ç

%%.0.1 Check bad channelsç
session = gui_session;
% if they have been updated, run
powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true);
powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true);
powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true);

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
ThetaChannels = []; % choose theta channels (ideally SLM)
SWChannels = []; % choose slow wave channels (ideally superficial cortex)
ignoretime = [];
% if EMG is not been quantified correctly, try discarting channels with
bz_EMGFromLFP(pwd,'rejectChannels',[],'overwrite', true,'ignoretime',ignoretime);

targetFile = (dir('*optogeneticPulses.events.mat'));
if ~isempty(targetFile)
    pulses = importdata(targetFile.name);
else
    pulses.stimulationEpochs = [];
end
ignoretime = [pulses.stimulationEpochs; ignoretime]; % [pulses.stimulationEpochs; ignoretime]
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',ignoretime, 'overwrite', true, 'ThetaChannels', ThetaChannels, 'SWChannels', SWChannels,'rejectChannels',[]);

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

ExcludeIntervals = [];
rippleChannel =29;
SWChannel = 6;
eventSpikeThreshold_shanks = [1 2 3 4 5]; % which shanks will be accounted for the spike threshold 
<<<<<<< HEAD
rippleMasterDetector_threshold = [1.5 3.5]; % [1.5 3.5]
eventSpikeThreshold = 1.5; % .5
=======
rippleMasterDetector_threshold = [.75 1.5]; % [1.5 3.5]
eventSpikeThreshold = 1.2; % .5
>>>>>>> e577f7ad5fa1dfc3a78c07fae9d009be16fdd81d
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'skipStimulationPeriods',false,'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold_shanks', eventSpikeThreshold_shanks,'eventSpikeThreshold',eventSpikeThreshold,'excludeIntervals',ExcludeIntervals); 
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
% psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);
getSpikesRank('events','ripples');

%% 5. Theta epochs
% Revise channel definition, theta band in thetaEpochs.png and cells
% rhytmicity. If bad, you can change useCSD_for_theta_detection to false,
% or change powerThreshold, even the channel

<<<<<<< HEAD
channel = 13;
=======
channel = 28;
>>>>>>> e577f7ad5fa1dfc3a78c07fae9d009be16fdd81d
useCSD_for_theta_detection = false;
powerThreshold = 1.6;% .8
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
% if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
%     file = dir([session.general.name,'.optogeneticPulses.events.mat']);
%     load(file.name);
% end
% excludeManipulationIntervals = optoPulses.stimulationEpochs;

cell_metrics = ProcessCellMetrics('session', session,'forceReload',true,'getWaveformsFromDat', true);

getACGPeak('force',true);

getAverageCCG('force',true);
    
getSpikesReturnPlot('force',true);

%% 9. Spatial modulation
% Open session folder with behaviour (with a video, normaly .avi, file). If
% no folder is present, you must run this code. You can play with the
% LED_threshold level (threshold for LED detection), and check if TTL
% channel are well defined. 
% NOS ESPERAMOS
LED_threshold = 0.98;
spikes = loadSpikes;
getSessionTracking('roiTracking','manual','forceReload',false,'LED_threshold',LED_threshold);
% only if the animal run a figure-eight maze behavior
%getSessionArmChoice('task','alternation','leftArmTtl_channel',2,'rightArmTtl_channel',3,'homeDelayTtl_channel',4);
behaviour = getSessionLinearize('forceReload',true);  
firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
firingTrialsMap = firingMapPerTrial('force',true);
spatialModulation = getSpatialModulation('force',true);

% if rReward or lReward dots are not clusterized in front of the IRSensor
% position, pray (work in progress)

%% 10. Checking analysis fiber
% %% 11. Spatial modulation
if ~any(ismember(excludeAnalysis, {'11',lower('spatialModulation')}))
    try
        spikes = loadSpikes;
        % getSessionTracking('roiTracking','manual','forceReload',false,'LED_threshold',LED_threshold,'convFact',tracking_pixel_cm,'leftTTL_reward',leftArmTtl_channel,'rightTTL_reward',rightArmTtl_channel);
        getSessionTracking('forceReload',false,'leftTTL_reward',leftArmTtl_channel,'rightTTL_reward',rightArmTtl_channel,'homeTtl',homeDelayTtl_channel,'tracking_ttl_channel',tracking_ttl_channel);
        try
            getSessionArmChoice('task','alternation','leftArmTtl_channel',leftArmTtl_channel,'rightArmTtl_channel',rightArmTtl_channel,'homeDelayTtl_channel',homeDelayTtl_channel,'use_manual_ttls',use_manual_ttls);
        catch
            warning('Performance in task was not computed! maybe linear maze?');
        end
        behaviour = getSessionLinearize('forceReload',false);  
        firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true,'speedThresh',0.1);
        placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
        firingTrialsMap = firingMapPerTrial('force',true,'saveMat',true);
        spatialModulation = getSpatialModulation('force',true);
    catch
        warning('Not possible to run spatial modulation...');
    end

    try 
        behaviour = getSessionLinearize;
        psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',100,'eventType','lReward','saveMat',false,...
            'minNumberOfPulses',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01,'raster_time',[-2 2]);
        psth_rReward = spikesPsth([behaviour.events.rReward],'numRep',100,'eventType','rReward','saveMat',false,...
            'minNumberOfPulses',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'raster_time',[-2 2]);
        psth_reward = spikesPsth([behaviour.events.lReward; behaviour.events.rReward],'numRep',100,'eventType','reward','saveMat',false,...
            'minNumberOfPulses',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'raster_time',[-2 2]);
        
        if all(isnan(behaviour.events.startPoint))
            behaviour.events.startPoint = NaN;
        end
        if all(isnan(behaviour.events.intersection))
            behaviour.events.intersection = NaN;
        end
        psth_intersection = spikesPsth([behaviour.events.intersection],'numRep',100,'eventType','intersection','saveMat',false,...
            'minNumberOfPulses',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'raster_time',[-2 2]);
        psth_startPoint = spikesPsth([behaviour.events.startPoint],'numRep',100,'eventType','startPoint','saveMat',false,...
            'minNumberOfPulses',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'raster_time',[-2 2]);

        behaviour.psth_lReward = psth_lReward;
        behaviour.psth_rReward = psth_rReward;
        behaviour.psth_reward = psth_reward;
        behaviour.psth_intersection = psth_intersection;
        behaviour.psth_startPoint = psth_startPoint; 
        behavior = behaviour; % british to american :)
        save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior','-v7.3');
    catch
        warning('Psth on behaviour events was not possible...');
    end

    % Fiber behaviour analysis
    try
        lReward_fiber = fiberPhotometryModulation_temp([behaviour.events.lReward],'eventType','lReward','saveMat',false);
        rReward_fiber = fiberPhotometryModulation_temp([behaviour.events.rReward],'eventType','rReward','saveMat',false);
        reward_fiber = fiberPhotometryModulation_temp([behaviour.events.lReward; behaviour.events.rReward],'eventType','reward','saveMat',false);

        intersection_fiber = fiberPhotometryModulation_temp([behaviour.events.intersection],'eventType','intersection','saveMat',false,'savePlotAs','intersection');
        startPoint_fiber = fiberPhotometryModulation_temp([behaviour.events.startPoint],'eventType','startPoint','saveMat',false,'savePlotAs','startPoint');

        fiber_behavior.lReward = lReward_fiber;
        fiber_behavior.rReward = rReward_fiber;
        fiber_behavior.reward = reward_fiber;

        fiber_behavior.intersection = intersection_fiber;
        fiber_behavior.startPoint = startPoint_fiber;

        save([basenameFromBasepath(pwd),'.behavior_fiber.events.mat'],'fiber_behavior');
    catch
        warning('No fiber recording in this session...');
    end

    try
        speedCorr = getSpeedCorr('numQuantiles',20,'force',true);
    catch
        warning('Speed\rate correlation analysis was not possible!');
    end
end


%% 10. Summary per cell
% Check how many cells have been included for the further analysis
% Also, if substancial changes were done before, run again!!
% (Summary_cell_1, 2 etc). If bad cells, run again. 
plotSummary('showTagCells',true);
close all; exit

% -------------------------------------------------------------------------



