
%%



digitalIn = getDigitalIn('force',true,'maze_in_virtual_channels',true);
optogeneticResponses = getOptogeneticResponse('force',true);

uLEDResponse = getuLEDResponse('force',true);
uLEDPulses = getuLEDPulses('force',true);



%
spikeFeatures('skipStimulationPeriods',false);
[phaseMod] = computePhaseModulation('skipStimulationPeriods',false);
session = loadSession;
cell_metrics = ProcessCellMetrics('session', session,'forceReload',true,'getWaveformsFromDat',true); % after CellExplorar
getAverageCCG('force',true,'skipStimulationPeriods',false);


TargetNeuroninfo = getSpikeTriggeredPulses;

