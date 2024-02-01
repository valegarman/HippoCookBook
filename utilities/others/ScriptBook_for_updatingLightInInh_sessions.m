
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
% 1. Copy ripple file (basename.ripples.events.mat) from NEURAL, and
% overwrite
% 2. Run this code;
targetFile = dir("*ripples.events.mat"); load(targetFile.name);
ripples.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel + 1;
ripples = computeRippleStats('ripples',ripples,'rippleChannel',ripples.detectorinfo.detectionchannel);
session = loadSession;
save([session.general.name , '.ripples.events.mat'],'ripples');
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
getSpikesRank('events','ripples');
[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel);
computeCofiringModulation;

close all;
[uLEDResponses] = getuLEDResponse('force',true);
% script_tempWin_LightInInh;
close all;
indexNewSession('copyFiles', true);

%% version for no pulses
pulses = getAnalogPulses('manualThr',true,'force',true); % 1-index   
close all;
getDigitalIn('force', true);
% uLEDPulses = getuLEDPulses('Current',3,'force',true,'ledLayout','ledLayoutScience2022');

% important, in the Extracellular tab, update Electrode groups and Spike groups witht following info:
% Group         Channel
% 1             25 21 24 28 26 22 23 27
% 2             29 17 20 32 30 18 19 31
% 3             16 4 1 13 15 3 2 14
% 4             12 8 5 9 11 7 6 10
processSession('digital_optogenetic_channels',[],'analog_optogenetic_channels',[],...
    'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false,...
    'bazler_ttl_channel',1,'leftArmTtl_channel',3,'leftArmTtl_channel',4,'useCSD_for_theta_detection',false);
close all
% for ripple detection, follow this steps:
% 1. Copy ripple file (basename.ripples.events.mat) from NEURAL, and
% overwrite
% 2. Run this code;
targetFile = dir("*ripples.events.mat"); load(targetFile.name);
ripples.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel + 1;
ripples = computeRippleStats('ripples',ripples,'rippleChannel',ripples.detectorinfo.detectionchannel);
session = loadSession;
save([session.general.name , '.ripples.events.mat'],'ripples');
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
getSpikesRank('events','ripples');
[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel);
computeCofiringModulation;

close all;
indexNewSession('copyFiles', true);
