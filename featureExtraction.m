function spikeFeatures = featureExtraction(varargin)
%
%       [spikeFeatures] = getSpikeFeatures(varargin)
%       
% Extract spike features that then will be used to separate neurons in
% different classes.
% 
% <OPTIONALS>
% basepath
% spikes                spikes struct  
% units ID              type of units to get features from (pyramidal,
%                       narror interneuron, wide interneuron).
% cell_metrics          cell_info struct. Options: 'pyr', 'int', nwint,
%                       wwint
% thetaModulation
% gammaModulation
% 
%
% OUTPUT
% spikeFeatures struct
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@bz_isCellInfo);
addParameter(p,'unitsID','all',@isstr);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'SW',[],@isstruct);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'SWpassband',[2 10], @isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
cell_metrics = p.Results.cell_metrics;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
unitsID = p.Results.unitsID;
SW = p.Results.SW;
SWChannel = p.Results.SWChannel;
SWpassband = p.Results.SWpassband;

%% Session Template
session = sessionTemplate(basepath,'showGUI',false);
%% Spikes
disp('Loading Spikes...')
spikes = loadSpikes;

%% Cell_metrics
if  ~isempty(dir([session.general.name,'.cell_metrics.cellinfo.mat']))
    disp('cell_metrics detected ! Loading file.')
    file = dir([session.general.name,'.cell_metrics.cellinfo.mat']);
    load(file.name)
else
    cell_metrics = ProcessCellMetrics('session', session);
end

if ischar(unitsID) && strcmpi(unitsID,'all')
    unitsID = spikes.UID;
elseif ischar(unitsID)
    for i = 1:length(cell_metrics.putativeCellType)
        if strcmpi(cell_metrics.putativeCellType{i},'Pyramidal Cell')
            celltype(i) = 1;
        elseif strcmpi(cell_metrics.putativeCellType{i},'Narrow Interneuron')
            celltype(i) = 2;
        elseif strcmpi(cell_metrics.putativeCellType{i},'Wide Interneuron')
            celltype(i) = 3;
        else
            celltype(i) = 4;
        end
    end
    
    switch unitsID
        case 'pyr'
            unitsID = find(cellType == 1);
        case 'int'
            unitsID = find(cellType == 2 | cellType == 3);
        case 'nwint'
            unitsID = find(cellType == 2);
        case 'wwint'
            unitsID = find(cellType == 3);
    end
end
    
%% Get optogenetic responses to know if a cell is statistically responding to stimulation
try
    if ~isempty(dir([session.general.name,'.optogeneticResponse.cellinfo.mat']))
        disp('Optogenetic Response detected. Loading file !');
        file = dir([session.general.name,'.optogeneticResponse.cellinfo.mat']);
        load(file.name);
    else
        optogeneticResponses = getOptogeneticResponse('numRep',100);
    end
catch
    disp('Not possible to get Optogenetic Responses..');
end
    
%% Phase Locking of Spikes to SharpWaves
% Loading SharpWaves
if isempty(SW)
    if ~isempty(dir([session.general.name,'.sharpwaves.events.mat']))
        disp('sharpwaves detected. Loading file...')
        file = dir([session.general.name,'.sharpwaves.events.mat']);
        load(file.name)
    else
        disp('sharpwaves not found. Running rippleMasterDetector');
        [ripples,SW] = rippleMasterDetector('SWChannel',SWChannel);
    end
end

if isempty(SWChannel) && ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
    file = dir([[session.general.name,'.hippocampalLayers.channelinfo.mat']]);
    load(file.name);
end
lfpSW = getLFP(SWChannel);
SWMod = SWphaseModulation(spikes,lfpSW,SWpassband,'intervals',SW.timestamps,'useThresh',false,'useMinWidth',false,'plotting',true);


%%
spikeFeatures = [];
spikeFeatures.ID = unitsID;
for i = 1:length(unitsID)
    % Spike events based metrics
    spikeFeatures.firingRate(i) = cell_metrics.firingRate(unitsID(i));
    % Waveform Based metrics
    spikeFeatures.troughtoPeak(i) = cell_metrics.troughToPeak(i);
    spikeFeatures.troughtoPeakDerivative(i) = cell_metrics.troughToPeakDerivative(i);
end




%% OUTPUT
spikeFeatures = [];
end

