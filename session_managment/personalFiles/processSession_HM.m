function [] = processSession_HM(varargin)

% [] = processSession_HM(varargin)

% This function runs standard analysis
% Based on index_session_script_InterneuronsLibrary and indexNewSession by MV 2020
%  1. 'sessionTemplate'         Runs sessionTemplate
%  2. 'loadSpikes'              Remove previous cellinfo.spikes.mat and computes spikes again (manual clustered)
%  3. 'cureAnalogPulses'        Analog pulses re-detection with threshold inspection
%  4. 'spikesFeatures'          Runs spike features: Light responses, if available, and ACG and waveform
%  5. 'checkSleep'              Check sleep score
%  6. 'powerProfiles'           Recompute power Profiles, considering bad channels now
%  7. 'eventsModulation'        Runs brain events modulation: i) Up and downs; ii) Ripples and iii) Theta intervals
%  8. 'phaseModulation'         Computes phase modulation for theta, gamma and ripples
%  9. 'cellMetrics'             Gets cell metrics (from Cell Explorer, and a few new)   
% 10. 'lfpAnalysis'             Computes power spectrum and coherence for
%                               different channels
% 11. 'summary'                 Makes cell/session sumary
%
% Note: to exclude any analysis, use 'excludeAnalysis' and provide the name
% or number of the section analysis to exlucde, example: processSession('excludeAnalysis',{'cureAnalogPulses', 'getHippocampalLayers', 8})
%
% Pablo Abad 2022, based on indexNewSession (now indeNewSession only indexes session on github)

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'project',[],@ischar); % SocialProject / SubiculumProject / MK801Project
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'rejectChannels',[],@isnumeric); % 0-index
addParameter(p,'force_analogPulsesDetection',true,@islogical);
addParameter(p,'force_loadingSpikes',true,@islogical);
addParameter(p,'excludeManipulationIntervals',[],@isnumeric);
addParameter(p,'rippleChannel',[],@isnumeric);% manually selecting ripple Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'SWChannel',[],@isnumeric); % manually selecting SW Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'thetaChannel',[],@isnumeric);
addParameter(p,'CA1Channel',[],@isnumeric); % manually selecting CA1 Channel / Can be also done by brain regions
addParameter(p,'SUBChannel',[],@isnumeric); % same as CA1 Channel
addParameter(p,'PFCChannel',[],@isnumeric); % same as CA1 Channel
addParameter(p,'digital_optogenetic_channels',[],@isnumeric);
addParameter(p,'analog_optogenetic_channels',[],@isnumeric);
addParameter(p,'promt_hippo_layers',false,@islogical);
addParameter(p,'manual_analog_pulses_threshold',false,@islogical);
addParameter(p,'bazler_ttl_channel',[],@isnumeric);
addParameter(p,'anymaze_ttl_channel',[],@isnumeric);
addParameter(p,'leftArmTtl_channel',2,@isnumeric)
addParameter(p,'rightArmTtl_channel',3,@isnumeric)
addParameter(p,'homeDelayTtl_channel',4,@isnumeric)
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'anyMaze',true,@islogical);
addParameter(p,'pixelsPerCm',2.5,@isnumeric);
addParameter(p,'excludeAnalysis',[]); % 
addParameter(p,'useCSD_for_theta_detection',false,@islogical);
addParameter(p,'profileType','hippocampus',@ischar); % options, 'hippocampus' and 'cortex'
addParameter(p,'showTetrodes',true,@islogical);
addParameter(p,'twoHalvesAnalysis',true,@islogical);
addParameter(p,'gridAnalysis',false,@islogical);
addParameter(p,'randomization',false,@islogical);
addParameter(p,'tint',true,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
project = p.Results.project;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;
rejectChannels = p.Results.rejectChannels;
force_analogPulsesDetection = p.Results.force_analogPulsesDetection;
force_loadingSpikes = p.Results.force_loadingSpikes;
excludeManipulationIntervals = p.Results.excludeManipulationIntervals;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thetaChannel = p.Results.thetaChannel;
CA1Channel = p.Results.CA1Channel;
SUBChannel = p.Results.SUBChannel;
PFCChannel = p.Results.PFCChannel;
digital_optogenetic_channels = p.Results.digital_optogenetic_channels;
analog_optogenetic_channels = p.Results.analog_optogenetic_channels;
promt_hippo_layers = p.Results.promt_hippo_layers;
manual_analog_pulses_threshold = p.Results.manual_analog_pulses_threshold;
bazler_ttl_channel = p.Results.bazler_ttl_channel;
anymaze_ttl_channel = p.Results.anymaze_ttl_channel;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
anyMaze = p.Results.anyMaze;
pixelsPerCm = p.Results.pixelsPerCm;
excludeAnalysis = p.Results.excludeAnalysis;
useCSD_for_theta_detection = p.Results.useCSD_for_theta_detection;
profileType = p.Results.profileType;
showTetrodes = p.Results.showTetrodes;
twoHalvesAnalysis = p.Results.twoHalvesAnalysis;
gridAnalysis = p.Results.gridAnalysis;
randomization = p.Results.randomization;
tint = p.Results.tint;

% Deal with inputs
prevPath = pwd;
cd(basepath);

for ii = 1:length(excludeAnalysis)
    excludeAnalysis{ii} = num2str(excludeAnalysis{ii});
end
if length(excludeAnalysis) == 0
    excludeAnalysis = num2str(excludeAnalysis);
end
excludeAnalysis = lower(excludeAnalysis);

<<<<<<< HEAD
getSessionStimulation();

try
    quadrants = getSessionQuadrants('force',true);
    
    psth_1 = spikesPsth([pulses.timestampsOn{1}'],'numRep',100,'saveMat',false,...
                'min_pulsesNumber',1,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 1],'binSize',0.01, 'win_Z',[-3 -1]);
    
    
    psth_quadrants_correct = spikesPsth([quadrants.ts(quadrants.choice == 1)],'numRep',100,'saveMat',false,...
                'min_pulsesNumber',1,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 1],'binSize',0.01, 'win_Z',[-3 -1]);
            
    psth_quadrants_incorrect = spikesPsth([quadrants.ts(quadrants.choice == 0)],'numRep',100,'saveMat',false,...
                'min_pulsesNumber',1,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 1],'binSize',0.01, 'win_Z',[-3 -1]);
catch
    warning('Not possible to run getQuadrants...');
end
=======
% try
%     quadrants = getSessionQuadrants('force',true);
%     
%     psth_quadrants_correct = spikesPsth([quadrants.ts(quadrants.choice == 1)],'numRep',100,'saveMat',false,...
%                 'min_pulsesNumber',1,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 1],'binSize',0.01, 'win_Z',[-3 -1]);
%             
%     psth_quadrants_incorrect = spikesPsth([quadrants.ts(quadrants.choice == 0)],'numRep',100,'saveMat',false,...
%                 'min_pulsesNumber',1,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-1 1],'binSize',0.01, 'win_Z',[-3 -1]);
% catch
%     warning('Not possible to run getQuadrants...');
% end
>>>>>>> 412eb6c149e8ffa0f85f46b3aaaf6d9572b8d322

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
        
        if ~isfield(session.analysisTags,'rippleChannel')
            session.analysisTags.rippleChannel = rippleChannel;
        end
        if ~isfield(session.analysisTags,'SWChannel')
            session.analysisTags.SWChannel = SWChannel;
        end
        if ~isfield(session.analysisTags,'thetaChannel')
            session.analysisTags.thetaChannel = thetaChannel;
        end
               
        save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
    catch
        warning('it seems that CellExplorer is not on your path');
    end

    session = sessionTemplate(basepath,'showGUI',true);

    selectProbe('force',true,'showTetrodes',showTetrodes); % choose probe
end

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
end

%% 3. Analog pulses detection
% if ~any(ismember(excludeAnalysis, {'3',lower('cureAnalogPulses')}))
%     if force_analogPulsesDetection || isempty(dir([session.general.name,'_original.dat']))
%         disp('Getting analog Pulses...')
%         pulses = getAnalogPulses('analogChannelsList',analog_optogenetic_channels,'manualThr',manual_analog_pulses_threshold,'overwrite',force_analogPulsesDetection); % 1-index
%     else
%         try
%             if ~isempty(dir([session.general.name,'.pulses.events.mat']))
%                 file = dir([session.general.name,'.pulses.events.mat']);
%                 load(file.name);
%             end
%         catch
%             warning('Problems with analogPulses');
%         end
%     end
% end

%% 4. Spike Features
% 4.1 Light responses, if available
if ~any(ismember(excludeAnalysis, {'4',lower('spikesFeatures')}))
%     optogeneticResponses = getOptogeneticResponse('numRep',500,'force',true);
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
%     TheStateEditor_temp(session.general.name);
    bz_ThetaStates(pwd);
end

%% 6. Power Profiles
if ~any(ismember(excludeAnalysis, {'6',lower('powerProfiles')}))
    powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'forceDetect',true);
    powerProfile_gamma = powerSpectrumProfile(gamma_bandpass,'showfig',true,'forceDetect',true);
    powerProfile_hfo = powerSpectrumProfile(hfo_bandpass,'showfig',true,'forceDetect',true);
end

%% 7. Check Brain Events
if ~any(ismember(excludeAnalysis, {'8',lower('eventsModulation')}))
    % Trying changes in detecUD_temp
    % 8.1 Up and downs
    UDStates = detectUD('plotOpt', true,'forceDetect',true','NREMInts','all');
    psthUD = spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true);

    % 8.2 Ripples
    ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'removeOptogeneticStimulation',true,'eventSpikeThreshold',false); % [1.5 3.5]
    psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);

    % 8.3 Theta intervals
    thetaEpochs = detectThetaEpochs('force',true,'useCSD',useCSD_for_theta_detection,'channel',thetaChannel);
end

%% 8. Phase Modulation
if ~any(ismember(excludeAnalysis, {'9',lower('phaseModulation')}))
    % LFP-spikes modulation
    [phaseMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'thetaChannel',thetaChannel,'lgammaChannel',lgammaChannel,'hgammaChannel',hgammaChannel);
    computeCofiringModulation;
end

%% 9. Cell metrics
% Exclude manipulation intervals for computing CellMetrics
if ~any(ismember(excludeAnalysis, {'10',lower('cellMetrics')}))
    if strcmpi(profileType,'hippocampus')
        session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]
    elseif strcmpi(profileType,'cortex')
        session = assignBrainRegion('showPowerProfile','hfo','showEvent','slowOscilations','eventTwin',[-.5 .5]); % hfo slowOscilations [-.5 .5]
    end

    try
        if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
            file = dir([session.general.name,'.optogeneticPulses.events.mat']);
            load(file.name);
        end
            excludeManipulationIntervals = optoPulses.stimulationEpochs;
    catch
        excludeManipulationIntervals = [];
        warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
    end
    cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'},'forceReload',true);

    getACGPeak('force',true);

    getAverageCCG('force',true);
    getAverageCCGPerSubSession('force',true);
    getSpikesReturnPlot('force',true);
%     computeAverageCCG('force',true);
    
    
end

%% 10. lfp Analysis
if ~any(ismember(excludeAnalysis,{'11',lower('lfpAnalysis')}))
    cohgram = computeCohgram('force',true);
end

imagination1 = computeImagination('channel',6);


semanticWords1 = computeSemanticWords('channel',1);
semanticWords2 = computeSemanticWords('channel',2);
semanticWords3 = computeSemanticWords('channel',3);
semanticWords4 = computeSemanticWords('channel',4);
semanticWords5 = computeSemanticWords('channel',5);

semanticWords9 = computeSemanticWords('channel',9);


onomatopeyas1 = computeOnomatopeyas('channel',2);
onomatopeyas1 = computeOnomatopeyas('channel',6);
onomatopeyas2 = computeOnomatopeyas('channel',7);
onomatopeyas3 = computeOnomatopeyas('channel',8);
onomatopeyas4 = computeOnomatopeyas('channel',10);
onomatopeyas5 = computeOnomatopeyas('channel',11);


file = dir('*pulses.events.mat');
load(file.name);

file = dir('*semantic_words.mat');
load(file.name);
indexes = cell2mat(words(:,2));

% PSTH
behavior = [];
psth_general = spikesPsth([pulses.timestampsOn{1}'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_visualObject = spikesPsth(pulses.timestampsOn{1}(indexes == 1)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_spatial = spikesPsth(pulses.timestampsOn{1}(indexes == 2)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_abstract = spikesPsth(pulses.timestampsOn{1}(indexes == 3)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_actions = spikesPsth(pulses.timestampsOn{1}(indexes == 4)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

% psth_familiar = spikesPsth(pulses.timestampsOn{9}(indexes == 5)','numRep',100,'saveMat',false,...
%     'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

behavior.psth_general = psth_general;
behavior.psth_visualObject = psth_visualObject;
behavior.psth_spatial = psth_spatial;
behavior.psth_abstract = psth_abstract;
behavior.psth_actions = psth_actions;
% behavior.psth_familiar = psth_familiar;

save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior');

% ONOMATOPEYAS
file = dir('*onomatopeyas_3.mat');
load(file.name);

file = dir('*pulses.events.mat');
load(file.name);
indexes = cell2mat(onomatopeyas(:,2:3));

onomat = [];

psth_primer = spikesPsth([pulses.timestampsOn{8}(indexes(:,1) == 1)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_stim = spikesPsth(pulses.timestampsOn{8}(indexes(:,1) == 2)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_coherent = spikesPsth([pulses.timestampsOn{8}(indexes(:,1) == 2 & indexes(:,2) == 1)'],'numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

psth_noncoherent = spikesPsth(pulses.timestampsOn{8}(indexes(:,1) == 2 & indexes(:,2) == 0)','numRep',100,'saveMat',false,...
    'min_PulsesNumber',5,'winSize',6,'binSize', 0.01,'event_ints',[0.8 1.5],'baseline_ints',[-0.4 -0.1],'winSizePlot',[-1 2]);

onomat.psth_primer = psth_primer;
onomat.psth_stim = psth_stim;
onomat.psth_coherent = psth_coherent;
onomat.psth_noncoherent = psth_noncoherent;

save([basenameFromBasepath(pwd) '.onomatopeyas.cellinfo.mat'],'onomat');

%% 11. Summary per cell
if ~any(ismember(excludeAnalysis, {'14',lower('summary')}))
    plotSummary_pablo();   
end

cd(prevPath);
end