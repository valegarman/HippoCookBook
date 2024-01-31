
%% ScriptBook_for_updatingLightInInh_sessions
% Marta Picco and Manu Valero, 2023
pulses = getAnalogPulses('manualThr',true,'force',true); % 1-index   
close all;
getDigitalIn('force', true);
uLEDPulses = getuLEDPulses('Current',3,'force',true,'ledLayout','ledLayoutScience2022');

% important, in the Extracellular tab, update Electrode groups and Spike groups witht following info:
% Group         Channel
% 1             25 21 24 28 26 22 23 27
% 2             29 17 20 32 30 18 19 31
% 3             16 4 1 13 15 3 2 14
% 4             12 8 5 9 11 7 6 10
processSession('digital_optogenetic_channels',[11 12 13 14 15 16],'analog_optogenetic_channels',[3 4 5 6 7 8],...
    'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false,...
    'bazler_ttl_channel',1,'leftArmTtl_channel',3,'leftArmTtl_channel',4,'useCSD_for_theta_detection',false);
close all
% for ripple detection, follow this steps:
% 1. Copy ripple file from NEURAL
% 2. Run this code;
targetFile = dir("*ripples.events.mat"); load(targetFile.name);
ripples.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel + 1;
ripples = computeRippleStats('ripples',ripples,'rippleChannel',ripples.detectorinfo.detectionchannel);
session = loadSession;
save([session.general.name , '.ripples.events.mat'],'ripples');
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
getSpikesRank('events','ripples');

[uLEDResponses] = getuLEDResponse('force',true);
script_tempWin_LightInInh;
indexNewSession('copyFiles', true);

%%  for sessions with manipulations (CNO and DMSO)
session = loadSession;
session = gui_session(session); % explore once last time! :)

uLEDResponses = getuLEDResponse('force',true);
uLEDResponses_post = getuLEDResponse('force',true,'restrict_to_manipulation',true);

getOptogeneticResponse('numRep',500,'force',true);
getOptogeneticResponse('numRep',500,'force',true,'restrict_to_manipulation',true);

powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true);
powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true);
powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true);
powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true,'minNumberOfPulses',10);
spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);

spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);

computePhaseModulation;
computePhaseModulation('restrict_to_manipulation',true);

if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
    load(file.name);
end
excludeManipulationIntervals = optoPulses.stimulationEpochs;
ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'restrictToIntervals',...
    uLEDResponses.restricted_interval,'manualAdjustMonoSyn',false);
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,...
    'restrictToIntervals',uLEDResponses.restricted_interval,'manualAdjustMonoSyn',false,'saveAs','cell_metrics_post');

getACGPeak('force',true);
getACGPeak('force',true,'restrict_to_manipulation',true);

getAverageCCG('force',true);
getAverageCCG('force',true,'restrict_to_manipulation',true);
    
getSpikesReturnPlot('force',true);
getSpikesReturnPlot('force',true,'restrict_to_manipulation',true);

getSpeedCorr('numQuantiles',20,'force',true);
getSpeedCorr('numQuantiles',20,'force',true, 'restrict_to_manipulation',true);

