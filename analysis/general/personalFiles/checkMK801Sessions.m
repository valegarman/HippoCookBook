%% ======== CHECKING MK801/GLUN3 SESSIONS ================

%% HPS22_100621_sess26 MK801 (Wildtype)

bpath = 'F:\data\HPS22\HPS22_100621_sess26';
cd(bpath);

session = loadSession();

SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[13257 str2num(session.general.duration)-1], 'overwrite', true);

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain region
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',[13257 str2num(session.general.duration)-1],'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[13257 str2num(session.general.duration)-1],'excludeMetrics',{'deepSuperficial'},'forceReload',true,'manualAdjustMonoSyn',false);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% HPS22_150621_sess27 (Wildtype)  VEHICLE (RERUN ANALYSIS)

bpath = 'F:\data\HPS22\HPS22_150621_sess27';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
close(gcf);plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');


%% HPS22_160621_sess28 (Wildtype) KETAMINE Animal dead when injected

bpath = 'F:\data\HPS22\HPS22_160621_sess28';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% HPS23_090621_sess9 (Wildtype) MK801

bpath = 'F:\data\HPS23\HPS23_090621_sess9';
cd(bpath);

session = loadSession(bpath);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);
            

%% HPS23_160621_sess11 (Wildtype) KETAMINE

bpath = 'F:\data\HPS23\HPS23_160621_sess11';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% HPS23_110621_sess10 (Wildtype) VEHICLE

bpath = 'F:\data\HPS23\HPS23_110621_sess10';
cd(bpath);

session = loadSession();

SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[], 'overwrite', true);

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true,'manualAdjustMonoSyn',false);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% HPS24_280621_sess5 (GLUN3) VEHICLE

bpath = 'F:\data\HPS24\HPS24_280621_sess5';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% HPS24_290621_sess6 (GLUN3) MK801

bpath = 'F:\data\HPS24\HPS24_290621_sess6';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% HPS24_230621_sess4 (GLUN3) KETAMINE

bpath = 'F:\data\HPS24\HPS24_230621_sess4';
cd(bpath);

session = loadSession();

SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[7257 7586], 'overwrite', true);

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% HPS25_050721_sess4 (GLUN3) MK801

bpath = 'F:\data\HPS25\HPS25_050721_sess4';
cd(bpath);

session = loadSession();

SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[814 874], 'overwrite', true);

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% HPS25_240621_sess2 (GLUN3) KETAMINE

bpath = 'F:\data\HPS25\HPS25_240621_sess2';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO430_111021_sess1 (GLUN3) VEHICLE

bpath = 'F:\data\IPO430\IPO430_111021_sess1';
cd(bpath);

session = loadSession();

SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',[4270 6922], 'overwrite', true);

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
file = dir(['*sharpwaves.events.mat']); load(file.name);

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',SW.detectorinfo.detectionchannel);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% IPO430_131021_sess2 (GLUN3) KETAMINE
bpath = 'F:\data\IPO430\IPO430_131021_sess2';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% IPO430_151021_sess3 (GLUN3) VEHICLE

bpath = 'F:\data\IPO430\IPO430_151021_sess3';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO429_111021_sess1 (Wild type) VEHICLE

bpath = 'F:\data\IPO429\IPO429_111021_sess1';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO429_131021_sess2 (Wild type) KETAMINE

bpath = 'F:\data\IPO429\IPO429_131021_sess2';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO429_151021_sess3 (Wild type) MK801

bpath = 'F:\data\IPO429\IPO429_151021_sess3';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});


% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO149_231021_sess3 (GLUN3) MK801

bpath = 'F:\data\IPO149\IPO149_231021_sess3';
cd(bpath);

session = loadSession();
% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO149_251021_sess4 (GLUN3) KETAMINE

bpath = 'F:\data\IPO149\IPO149_251021_sess4';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO149_281021_sess5 (GLUN3) VEHICLE

bpath = 'F:\data\IPO149\IPO149_281021_sess5';
cd(bpath);

session = loadSession();

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO150_261021_sess2 (GLUN3) MK801
bpath = 'F:\data\IPO150\IPO150_091121_sess4';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);

%% IPO150_091121_sess4 (GLUN3) VEHICLE

bpath = 'F:\data\IPO150\IPO150_091121_sess4';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% IPO150_201921_sess3 KETAMINE

bpath = 'F:\data\IPO150\IPO150_291021_sess3';
cd(bpath);

session = loadSession();

% Checking States
TheStateEditor(session.general.name);

% Rerun thetaEpochs and modulation
% 8.3 Theta intervals
file = dir(['*thetaEpochs.states.mat']); load(file.name);
thetaEpochs = detectThetaEpochs('force',true,'useCSD',false,'channel',thetaEpochs.channel);

% LFP-spikes modulation
file = dir(['*ripples.events.mat']); load(file.name);
try
    file = dir(['*sharpwaves.events.mat']); load(file.name);
catch
end

[phaseMod] = computePhaseModulation('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',13);
computeCofiringModulation;

 % LFP-spikes modulation per subsession
[phaseModSubSession] = computePhaseModulationPerSubSession('rippleChannel',ripples.detectorinfo.detectionchannel,'SWChannel',[]);

% Brain Regions
session = assignBrainRegion('showPowerProfile','theta','showEvent','ripples','eventTwin',[-0.05 0.05]); % hfo slowOscilations [-.5 .5]

% Rerun ProcessCellmetrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});

% Cell Metrics per subsession
cell_metrics_SubSession = ProcessCellMetricsPerSubSession('session', session,'excludeIntervals',[],'excludeMetrics',{'deepSuperficial'},'forceReload',true);
    
% Checking CellExplorer
cell_metrics = loadCellMetrics('basepath',pwd);
cell_metrics = CellExplorer('metrics',cell_metrics);

% Plotting again
plotSummary_pablo('excludePlot',{'spatialModulation'},'saveAs','Checked');
close(gcf);


%% IPO447_181121_sess1 VEHICLE


%% IPO447_241121_sess4 KETAMINE


