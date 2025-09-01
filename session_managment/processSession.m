function [] = processSession(varargin)

% [] = indexSession_InterneuronsLibrary(varargin)

% This function runs standard analysis
% Based on index_session_script_InterneuronsLibrary and indexNewSession by MV 2020
%  1. 'sessionTemplate'         Runs sessionTemplate
%  2. 'loadSpikes'              Remove previous cellinfo.spikes.mat and computes spikes again (manual clustered)
%  3. 'cureAnalogPulses'        Analog pulses re-detection with threshold inspection
%  4. 'spikesFeatures'          Runs spike features: Light responses, if available, and ACG and waveform
%  5. 'checkSleep'              Check sleep score
%  6. 'powerProfiles'           Recompute power Profiles, considering bad channels now
%  7. 'getHippocampalLayers'    Define semi-automatically hippocampal Layers
%  8. 'eventsModulation'        Runs brain events modulation: i) Up and downs; ii) Ripples and iii) Theta intervals
%  9. 'phaseModulation'         Computes phase modulation for theta, gamma and ripples
% 10. 'cellMetrics'             Gets cell metrics (from Cell Explorer, and a few new)                       
% 11. 'spatialModulation'       Process spatial modulation analysis, behavioural events and speed
% 12. 'summary'                 Makes cell/session sumary
%
% Note: to exclude any analysis, use 'excludeAnalysis' and provide the name
% or number of the section analysis to exlucde, example: processSession('excludeAnalysis',{'cureAnalogPulses', 'getHippocampalLayers', 8})
%
% MV 2022, based on indexNewSession (now indeNewSession only indexes session on github)

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'theta_bandpass',[5 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'rejectChannels',[],@isnumeric); % 0-index
addParameter(p,'force_analogPulsesDetection',true,@islogical);
addParameter(p,'force_loadingSpikes',true,@islogical);
addParameter(p,'excludePulsesIntervals',[],@isnumeric);
addParameter(p,'thetaChannel',[],@isnumeric);
addParameter(p,'theta_epochs_channel',[]);
addParameter(p,'rippleChannel',[],@isnumeric);% manually selecting ripple Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'SWChannel',[],@isnumeric); % manually selecting SW Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'digital_optogenetic_channels',[],@isnumeric);
addParameter(p,'analog_optogenetic_channels',[],@isnumeric);
addParameter(p,'promt_hippo_layers',true,@islogical);
addParameter(p,'manual_analog_pulses_threshold',false,@islogical);
addParameter(p,'bazler_ttl_channel',[],@isnumeric);
addParameter(p,'leftArmTtl_channel',[],@isnumeric);
addParameter(p,'rightArmTtl_channel',[],@isnumeric);
addParameter(p,'homeDelayTtl_channel',[],@isnumeric);
addParameter(p,'dlc_ttl_channel',[],@isnumeric);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'excludeAnalysis',[]); % 
addParameter(p,'profileType','hippocampus',@ischar); % options, 'hippocampus' and 'cortex'
addParameter(p,'rippleMasterDetector_threshold',[1.5 3.5],@isnumeric); % [1.5 3.5]
addParameter(p,'LED_threshold',0.98,@isnumeric);
addParameter(p,'createLegacySummaryFolder',false,@islogical);
addParameter(p,'useCSD_for_theta_detection',false,@islogical);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'selectProbe_automatic',false,@islogical);
addParameter(p,'use_manual_ttls',false,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;
rejectChannels = p.Results.rejectChannels;
force_analogPulsesDetection = p.Results.force_analogPulsesDetection;
force_loadingSpikes = p.Results.force_loadingSpikes;
excludePulsesIntervals = p.Results.excludePulsesIntervals;
thetaChannel = p.Results.thetaChannel;
theta_epochs_channel = p.Results.theta_epochs_channel;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
digital_optogenetic_channels = p.Results.digital_optogenetic_channels;
analog_optogenetic_channels = p.Results.analog_optogenetic_channels;
promt_hippo_layers = p.Results.promt_hippo_layers;
manual_analog_pulses_threshold = p.Results.manual_analog_pulses_threshold;
bazler_ttl_channel = p.Results.bazler_ttl_channel;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
dlc_ttl_channel = p.Results.dlc_ttl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
excludeAnalysis = p.Results.excludeAnalysis;
profileType = p.Results.profileType;
rippleMasterDetector_threshold = p.Results.rippleMasterDetector_threshold;
LED_threshold = p.Results.LED_threshold;
useCSD_for_theta_detection = p.Results.useCSD_for_theta_detection;
createLegacySummaryFolder = p.Results.createLegacySummaryFolder;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
selectProbe_automatic = p.Results.selectProbe_automatic;
useCSD_for_theta_detection = p.Results.useCSD_for_theta_detection;
use_manual_ttls = p.Results.use_manual_ttls;

% Deal with inputs
prevPath = pwd;
cd(basepath);

mkdir('SummaryFigures')
if createLegacySummaryFolder
    if exist('SummaryFigures') == 7
        d = strrep(strrep(string(datetime),' ','_'),':','_');
        cd('SummaryFigures\');
        legacyFolderName = strcat('legacy',strrep(strrep(string(datetime),' ','_'),':','_'),'_SummaryFigures');
        mkdir(legacyFolderName)
        allFigures = dir('*.png');
        for ii = 1:length(allFigures)
            movefile(allFigures(ii).name,[char(legacyFolderName) filesep allFigures(ii).name]);
        end
        cd(basepath);
    end
end

for ii = 1:length(excludeAnalysis)
    excludeAnalysis{ii} = num2str(excludeAnalysis{ii});
end
if length(excludeAnalysis) == 0
    excludeAnalysis = num2str(excludeAnalysis);
end
excludeAnalysis = lower(excludeAnalysis);

%% 1. Runs sessionTemplate
if ~any(ismember(excludeAnalysis, {'1',lower('sessionTemplate')}))
        try
            session = loadSession(basepath);
            session.channels = 1:session.extracellular.nChannels;
    
            if ~isfield(session, 'analysisTags')
                session.analysisTags = [];
            end
            if ~isfield(session.analysisTags,'digital_optogenetic_channels')
                session.analysisTags.digital_optogenetic_channels = digital_optogenetic_channels;
            end
            if ~isfield(session.analysisTags,'analog_optogenetic_channels')
                session.analysisTags.analog_optogenetic_channels = analog_optogenetic_channels;
            end
    
            if isempty(rejectChannels)
                rejectChannels = session.channelTags.Bad.channels; % 1-index
            end
            if ~isfield(session.analysisTags,'bazler_ttl_channel')
                session.analysisTags.bazler_ttl_channel = bazler_ttl_channel;
            end
            if ~isfield(session.analysisTags,'dlc_ttl_channel')
                session.analysisTags.dlc_ttl_channel = dlc_ttl_channel;
            end
            
            if ~isfield(session.analysisTags,'leftArmTtl_channel')
                session.analysisTags.leftArmTtl_channel = leftArmTtl_channel;
            end
            if ~isfield(session.analysisTags,'rightArmTtl_channel')
                session.analysisTags.rightArmTtl_channel = rightArmTtl_channel;
            end
            if ~isfield(session.analysisTags,'homeDelayTtl_channel')
                session.analysisTags.homeDelayTtl_channel = homeDelayTtl_channel;
            end
            
            if ~isfield(session.analysisTags,'homeDelayTtl_channel')
            end
            
            save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
        catch
            warning('it seems that CellExplorer is not on your path');
            
            session = sessionTemplate(basepath,'showGUI',true);
        end


    session = gui_session(session);

    selectProbe('force',true,'automatic', selectProbe_automatic); % choose probe
    close all
end

leftArmTtl_channel = session.analysisTags.leftArmTtl_channel;
rightArmTtl_channel = session.analysisTags.rightArmTtl_channel;
homeDelayTtl_channel = session.analysisTags.homeDelayTtl_channel;
tracking_ttl_channel = session.analysisTags.bazler_ttl_channel;

ints = [];
if restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif restrict_to_baseline
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [0 session.epochs{ii}.startTime];
            warning('Epoch with manipulations found! Restricting analysis to baseline interval!');
        end
    end
    if isempty(ints)
        ints = [0 Inf];
    end
else
    ints = [0 Inf];
end

restrict_ints = IntersectIntervals([ints; restrict_to]);


%% 2. Remove previous cellinfo.spikes.mat and computes spikes again (manual clustered)
if ~any(ismember(excludeAnalysis, {'2',lower('loadSpikes')}))
    disp('Loading Spikes...')
    session = loadSession;
    if force_loadingSpikes
        if ~isempty(dir([session.general.name ,'.spikes.cellinfo.mat']))
            file = dir([session.general.name ,'.spikes.cellinfo.mat']);
            delete(file.name);
        end
    end
    spikes = loadSpikes('forceReload',force_loadingSpikes);
    
    session = loadSession(basepath);
end

%% 3. Analog pulses detection
if ~any(ismember(excludeAnalysis, {'3',lower('cureAnalogPulses')}))
    if force_analogPulsesDetection || isempty(dir([session.general.name,'_original.dat']))
        disp('Getting analog Pulses...')
        pulses = getAnalogPulses('analogChannelsList',analog_optogenetic_channels,'manualThr',manual_analog_pulses_threshold,'force',force_analogPulsesDetection); % 1-index
    else
        try
            if ~isempty(dir([session.general.name,'.pulses.events.mat']))
                file = dir([session.general.name,'.pulses.events.mat']);
                load(file.name);
            end
        catch
            warning('Problems with analogPulses');
        end
    end
end

%% 4. Spike Features, and optogenetic responses
% 4.1 Light responses, if available
if ~any(ismember(excludeAnalysis, {'4',lower('spikesFeatures')}))
    % getOptogeneticResponse('numRep',500,'force',true,'restrict_to', restrict_ints,'digitalChannelsList',digital_optogenetic_channels);
    % 4.2 ACG and waveform
    spikeFeatures;
end

%% 5. Check Sleep Score
if ~any(ismember(excludeAnalysis, {'5',lower('checkSleep')}))
    targetFile = (dir('*optogeneticPulses.events.mat'));
    if ~isempty(targetFile)
        pulses = importdata(targetFile.name);
    else
        pulses.stimulationEpochs = [];
    end
    
    SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.stimulationEpochs, 'overwrite', true);
%     TheStateEditor(session.general.name);
    bz_ThetaStates(pwd);
end

%% 6. Power Profiles
if ~any(ismember(excludeAnalysis, {'6',lower('powerProfiles')}))
    powerSpectrumProfile(theta_bandpass,'showfig',true,'forceDetect',true,'restrict_to',restrict_ints);
    powerSpectrumProfile(gamma_bandpass,'showfig',true,'forceDetect',true,'restrict_to',restrict_ints);
    powerSpectrumProfile(hfo_bandpass,'showfig',true,'forceDetect',true,'restrict_to',restrict_ints);
end

%% 7. Getting Hippocampal Layers
if ~any(ismember(excludeAnalysis, {'7',lower('getHippocampalLayers')}))
    [hippocampalLayers] = getHippocampalLayers('force',true,'promt',promt_hippo_layers,'removeRipplesStimulation', false);
end

%% 8. Check Brain Events
if ~any(ismember(excludeAnalysis, {'8',lower('eventsModulation')}))
    % 8.1 Up and downs
    UDStates = detectUpsDowns('plotOpt', true,'forceDetect',true,'NREMInts','all','useparfor',true);
    psthUD = spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to',restrict_ints);
    getSpikesRank('events','upstates');

    % 8.2 Ripples

    ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'skipStimulationPeriods',true,'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold', false); 
    psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
    % psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);
    getSpikesRank('events','ripples');

    % 8.4 Fiber ripple analysis
    try
        ripples_fiber = fiberPhotometryModulation_temp([],'eventType','ripples');
    catch
        warning('No fiber recording in this session...');
    end

    % 8.3 Theta intervals
    cd(basepath);
    thetaEpochs = detectThetaEpochs('force',true,'useCSD',useCSD_for_theta_detection,'channel',theta_epochs_channel);

end


%% 9. Phase Modulation
if ~any(ismember(excludeAnalysis, {'9',lower('phaseModulation')}))
    % LFP-spikes modulation
    [phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',thetaChannel,'hgammaChannel',thetaChannel,'restrict_to',restrict_ints);
    % [phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'restrict_to',restrict_ints);
    computeCofiringModulation;
end

%% 10. Cell metrics
% Exclude manipulation intervals for computing CellMetrics
if ~any(ismember(excludeAnalysis, {'10',lower('cellMetrics')}))
    if strcmpi(profileType,'hippocampus')
        session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]);
    elseif strcmpi(profileType,'cortex')
        session = assignBrainRegion('showPowerProfile','hfo','showEvent','slowOscilations','eventTwin',[-.5 .5]); % hfo slowOscilations [-.5 .5]
    end
    
    if isempty(excludePulsesIntervals)
        try
            if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
                file = dir([session.general.name,'.optogeneticPulses.events.mat']);
                load(file.name);
            end
                excludePulsesIntervals = optoPulses.stimulationEpochs;
        catch
            warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
        end
    end

    session = loadSession;
        
    cell_metrics = ProcessCellMetrics('session', session,'forceReload',true,'excludeIntervals',excludePulsesIntervals,'excludeMetrics',{'deepSuperficial'}); % after CellExplorar
        
    getACGPeak('force',true);

    getAverageCCG('force',true);
    
    getSpikesReturnPlot('force',true);
end

%% 11. Spatial modulation
if ~any(ismember(excludeAnalysis, {'11',lower('spatialModulation')}))
    try
        spikes = loadSpikes;
        %getSessionTracking('roiTracking','manual','forceReload',false,'LED_threshold',LED_threshold,'convFact',tracking_pixel_cm,'leftTTL_reward',leftArmTtl_channel,'rightTTL_reward',rightArmTtl_channel);
        getSessionTracking('forceReload',true,'leftTTL_reward',leftArmTtl_channel,'rightTTL_reward',rightArmTtl_channel,'homeTtl',homeDelayTtl_channel,'dlc_ttl_channel',dlc_ttl_channel);
        try
            getSessionArmChoice('task','alternation','leftArmTtl_channel',leftArmTtl_channel,'rightArmTtl_channel',rightArmTtl_channel,'homeDelayTtl_channel',homeDelayTtl_channel,'use_manual_ttls',use_manual_ttls);
        catch
            warning('Performance in task was not computed! maybe linear maze?');
        end
        behaviour = getSessionLinearize('forceReload',true,'leftTtl',leftArmTtl_channel,'rightTtl',rightArmTtl_channel);  
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

        intersection_fiber = fiberPhotometryModulation_temp([behaviour.events.intersection],'eventType','intersection','saveMat',false);
        startPoint_fiber = fiberPhotometryModulation_temp([behaviour.events.startPoint],'eventType','startPoint','saveMat',false);

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

%% 12. ULED analysis 
try
    getuLEDPulses;
    getuLEDResponse;         
catch
    warning('Not possible to run getULEDResponse...');
end

%% 13. Summary per cell
if ~any(ismember(excludeAnalysis, {'12',lower('summary')}))
    plotSummary('showTagCells',false,'use_deltaThetaEpochs',false);
    
end

cd(prevPath);
end




