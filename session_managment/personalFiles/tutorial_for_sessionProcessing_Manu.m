%% Tutorial for session processing

% 1% First, transfer files from recording computer to analysis unit (wheter
%   a computer, NASS or Cloud Share site) by
%   'updateExpFolder({recordingPC_1, recordingPC_2, etc}, 'analysis unit')',
%   Example:

updateExpFolder({'V:\data\fCck1', 'Y:\fCck1'},'E:\data\fCck1');

% 2% Then, preprocess session (includes artifacts removal, median signal
%   removal, LFP and Kilosort, and running computeSessionSummary by 'batch_preprocessSession('basepath','sessionBasepath').
%   Example:
batch_preprocessSession('basepath','E:\wt3');

% <OPTIONAL> If summary was not processed, it can be run in batch by 'batch_preprocessSession'
batch_sessionSummary('basepath','G:\data\fPv4','cleanArtifacts',({65,[]}),'analogChannelsList',65,'digitalChannelsList',0);

% <OPTIONAL> If only a single session should be processed
preprocessSession('basepath',pwd,'analogChannelsList',[1:3],'spikeSort',true,'getPos',false, 'cleanArtifacts', true,...
                    'medianSubstr',true,'tracking_pixel_cm',NaN,'sessionSummary',true,'digitalChannelsList',[],'bazler_ttl_channel',[],'skipStimulationPeriods',false);

% 3% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[0],'analog_optogenetic_channels',[1],'promt_hippo_layers',true,'profileType','hippocampus');

% 5% Index session
indexNewSession('copyFiles',true);

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\valeg\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',false);

% run getOptogeneticResponses with uLEDs
pulses = getAnalogPulses('manualThr',true,'force',true); % 1-index
getDigitalIn;
uledPulses = getuLEDPulses('current',6,'force', true);
optogeneticResponses = getOptogeneticResponse('numRep',50,'force',true,...
    'analogChannelsList',[3 4 5 6 7 8],'digitalChannelsList',[11 12 13 14 15 16]);
getuLEDResponse('force',true);

%% code for humans
preprocessSession('basepath',pwd,'spikeSort',true,'getPos',false, ...
                    'medianSubstr',true);
% CLEAN SESSIONS MANUALLY BY PHY
processSession('promt_hippo_layers',true,'profileType','hippocampus','selectProbe_automatic', false);
cellTypes = cellTypeClassifier('modelType','fobrebrain5','score_cut_off',0,'imposeCellExplorerPyr',false);
indexNewSession('copyFiles',true);

offset = 37.5;
for ii = 1:10
    disp(offset) 
    offset = offset + 25;
end

linspace()