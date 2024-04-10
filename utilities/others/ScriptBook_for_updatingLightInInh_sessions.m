
%% ScriptBook_for_updatingLightInInh_sessions
% Marta Picco and Manu Valero, 2023

%% 1. Session with pulses 

% %%% getAnalogPulses generates for each channel the recording in order to
% % dselect the best channel to be used for thanalysis and save different
% % features about the pulses.
% 
% pulses = getAnalogPulses('manualThr',true,'force',true); % 1-index   
% getDigitalIn('force', true);
% uLEDPulses = getuLEDPulses('Current',3,'force',true,'ledLayout','ledLayoutScience2022');
% close all;
% 
% % ProcessSession is a big funtion that runs different analysis. Go inside
% % the function to check the different section.
% % important, in the Extracellular tab, update Electrode groups and Spike groups witht following info:
% % Group         Channel
% % 1             25 21 24 28 26 22 23 27
% % 2             29 17 20 32 30 18 19 31
% % 3             16 4 1 13 15 3 2 14
% % 4             12 8 5 9 11 7 6 10
% 
% processSession('digital_optogenetic_channels',[11 12 13 14 15 16],'analog_optogenetic_channels',[3 4 5 6 7 8],...
%     'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false,...
%     'bazler_ttl_channel',1,'leftArmTtl_channel',3,'leftArmTtl_channel',4,'useCSD_for_theta_detection',false);
% close all%%
% For ripple detection, follow this steps:
% % 1. Copy ripple file (basename.ripples.events.mat) from NEURAL, and
% % overwrite in the folder of the Session in your local memory
% % 2. Run this code;
% 
% targetFile = dir("*ripples.events.mat"); load(targetFile.name);
% ripples.detectorinfo.detectionchannel = ripples.detectorinfo.detectionchannel + 1;
% ripples = computeRippleStats('ripples',ripples,'rippleChannel',ripples.detectorinfo.detectionchannel);
% session = loadSession;
% save([session.general.name , '.ripples.events.mat'],'ripples');
% psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
% getSpikesRank('events','ripples');
% [phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel);
% computeCofiringModulation;
% close all;


[uLEDResponses] = getuLEDResponse('force',true); %questa funzione ha senso runnarla solo se si ha usato il LED
% script_tempWin_LightInInh;
close all;

indexNewSession('copyFiles', true);

%% 2.Session with no pulses 

pulses = getAnalogPulses('manualThr',true,'force',true); % 1-index   
close all;
getDigitalIn('force', true);


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

% For ripple detection, follow this steps:
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

%% 3. Sessions with manipulations (CNO and DMSO) after have already runned the section 1 or 2

% New window will open in order to express General, Epochs, Animal Subject,
% Extracellular, Tags, etc..

session = loadSession;
session = gui_session(session); % explore once last time! :)

%If no pulses, go directly to the Power spectrum Profuile function

% getuLEDResponse generates two plots showing, for each neuron detected, 
% how much and by which mLED is activated.The same result is shown using z_score

uLEDResponses = getuLEDResponse('force',true);
uLEDResponses_post = getuLEDResponse('force',true,'restrict_to_manipulation',true);

disp ('uLEDResponses computed!');

% getOptogeneticResponse generates N plots, as N the number of the mLED
% used in the experiment.For each mLED is showed the behaviour of all
% neurons detected trough a raster plot. 

getOptogeneticResponse('numRep',500,'force',true);
getOptogeneticResponse('numRep',500,'force',true,'restrict_to_manipulation',true);

disp ('OptogeneticResponse computed!');

% Three different frequency bands of interest: theta power ( 6-12 Hz), gamma power (20-100 Hz), RHF (100-500 Hz).
% powerSpectrumProfile, for each frequency range, generates a plot. 
% The graph shows for each shank, at the level of each channel present on that shank (y_axis),
% the power of the frequency band being analyzed (x_axis).  

powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true);
powerSpectrumProfile([6 12],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true);
powerSpectrumProfile([20 100],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true);
powerSpectrumProfile([100 500],'showfig',true,'forceDetect',true,'restrict_to_manipulation',true);

disp ('powerSpectrumProfile computed!');

% spikesPsth is a function that generates a Peri Stimolous Time Histogram
% around an event expressed in 'eventType' and in a time interval. 

spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true,'minNumberOfPulses',10);
spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);

spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10,'restrict_to_manipulation',true);

disp('Peri Stimolous Time Histograms computed!');

% computePhaseModulation compute the phase distribution of the firing of
% the neurons for difference frequencies bands. 

computePhaseModulation;
computePhaseModulation('restrict_to_manipulation',true);

disp('PhaseModulation computed!');
 
if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
    file = dir([session.general.name,'.optogeneticPulses.events.mat']);
    load(file.name);
end

% Getting waveforms of the current session and using the feature explained in CellMetrics website.
% ProcessCellMetrics generates plot about CellMetrics session summary 

% if no pulses, put this variable related to the pulses or to the uLED response equal to []

excludeManipulationIntervals = optoPulses.stimulationEpochs; 

ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,'restrictToIntervals',...
    uLEDResponses.restricted_interval,'manualAdjustMonoSyn',false);
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'forceReload',true,...
    'restrictToIntervals',uLEDResponses.restricted_interval,'manualAdjustMonoSyn',false,'saveAs','cell_metrics_post');

disp('CellMetrics computed!');

% getACGPeak compute an auto correlation of each neuron in log scale. To
% expand the part of interest. 

getACGPeak('force',true);
getACGPeak('force',true,'restrict_to_manipulation',true);


% getAverageCCG compute, for each neurons detected a cross correlation between itself and 
% all other neurons. The plot show how each neurons is efficiently involved in the network
% activity. 

getAverageCCG('force',true);
getAverageCCG('force',true,'restrict_to_manipulation',true); 

disp('ACGpeak and AvarageCCG computed!');

% getSpikesReturnPlot compute the time distance between one neuron to all
% other nurons using this information as coordinates in the cartesian
% plane.

getSpikesReturnPlot('force',true);
getSpikesReturnPlot('force',true,'restrict_to_manipulation',true);

disp('SpikesReturnPlot computed!');

% getSpeedCorr compute correlation between the speed of the animal and the firing rate
% of each neurons.

getSpeedCorr('numQuantiles',20,'force',true);
getSpeedCorr('numQuantiles',20,'force',true, 'restrict_to_manipulation',true);

disp ('Analysis completed!');
