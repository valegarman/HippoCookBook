function cell_metrics = ProcessCellMetricsPerSubSession(varargin)
%   This function calculates cell metrics for a given recording/session
%   Most metrics are single value per cell, either numeric or string type, but
%   certain metrics are vectors like the autocorrelograms or cell with double content like waveforms.
%   The metrics are based on a number of features: spikes, waveforms, PCA features,
%   the ACG and CCGs, LFP, theta, ripples and so fourth
%
%   Check the website of CellExplorer for more details: https://cellexplorer.org/
%
%   % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   INPUTS
%   % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   varargin (Variable-length input argument list; see below)
%
%   - Parameters defining the session to process - 
%   basepath               - 1. Path to session (base directory)
%   sessionName            - 3. Database sessionName
%   sessionID              - 4. Database numeric id
%   session                - 5. Session struct. Must contain a basepath
%
%   - Parameters for the processing - parameters.*
%   showGUI                - Show GUI dialog to adjust settings/parameters
%   metrics                - Metrics that will be calculated. A cell with strings
%                            Examples: 'waveform_metrics','PCA_features','acg_metrics','deepSuperficial',
%                            'monoSynaptic_connections','theta_metrics','spatial_metrics',
%                            'event_metrics','manipulation_metrics', 'state_metrics','psth_metrics'
%                            Default: 'all'
%   excludeMetrics         - Metrics to exclude. Default: 'none'
%   removeMetrics          - Metrics to remove (supports only deepSuperficial at this point)
%   keepCellClassification - logical. Keep existing cell type classifications
%   includeInhibitoryConnections - logical. Determines if inhibitory connections are included in the detection of synaptic connections
%   manualAdjustMonoSyn    - logical. Manually validate monosynaptic connections in the pipeline (requires user input)
%   restrictToIntervals    - time intervals to restrict the analysis to (in seconds)
%   excludeIntervals       - time intervals to exclude (in seconds)
%   excludeManipulationIntervals - logical. Exclude time intervals around manipulations (loads *.manipulation.mat files and excludes defined manipulation intervals)
%   ignoreEventTypes       - exclude .events files of specific types
%   ignoreManipulationTypes- exclude .manipulations files of specific types
%   ignoreStateTypes       - exclude .states files of specific types
%   showGUI                - logical. Show a GUI that allows you to adjust the input parameters/settings
%   forceReload            - logical. Recalculate existing metrics
%   forceReloadSpikes      - logical. Reloads spikes and other cellinfo structs
%   submitToDatabase       - logical. Submit cell metrics to database
%   saveMat                - logical. Save metrics to cell_metrics.mat
%   saveAs                 - name of .mat file
%   saveBackup             - logical. Whether a backup file should be created
%   summaryFigures         - logical. Plot summary figures
%   debugMode              - logical. Activate a debug mode avoiding try/catch 
%   transferFilesFromClusterpath - logical. Moves previosly generated files from clusteringpath to basepath (new file structure)
%   showFigures            - logical. if false, turns off the default plotting of different stages of processing 
%
% - Example calls:
%   cell_metrics = ProcessCellMetrics                             % Load from current path, assumed to be a basepath
%   cell_metrics = ProcessCellMetrics('session',session)          % Load session from session struct
%   cell_metrics = ProcessCellMetrics('basepath',basepath)        % Load from basepath
%   cell_metrics = ProcessCellMetrics('basepath',basepath,'includeInhibitoryConnections',true) % Load from basepath
%   cell_metrics = ProcessCellMetrics('basepath',basepath,'metrics',{'waveform_metrics'}) % Load from basepath
%
%   cell_metrics = ProcessCellMetrics('sessionName','rec1')       % Load session from database session name
%   cell_metrics = ProcessCellMetrics('sessionID',10985)          % Load session from database session id
%  
%   cell_metrics = ProcessCellMetrics('session', session,'fileFormat','nwb','showGUI',true); % saves cell_metrics to a nwb file
%   
%   % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   OUTPUT
%   % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   cell_metrics : structure described in details at: https://cellexplorer.org/datastructure/standard-cell-metrics/

%   By Peter Petersen
%   Last edited: 27-02-2021

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Parsing parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

p = inputParser;
addParameter(p,'sessionID',[],@isnumeric);
addParameter(p,'sessionName',[],@isstr);
addParameter(p,'session',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'metrics','all',@iscellstr);
addParameter(p,'excludeMetrics',{'none'},@iscellstr);
addParameter(p,'removeMetrics',{'none'},@isstr);
addParameter(p,'restrictToIntervals',[],@isnumeric);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'ignoreEventTypes',{'MergePoints'},@iscell);
addParameter(p,'ignoreManipulationTypes',{'cooling'},@iscell);
addParameter(p,'ignoreStateTypes',{'StateToIgnore'},@iscell);
addParameter(p,'excludeManipulationIntervals',true,@islogical);
addParameter(p,'metricsToExcludeManipulationIntervals',{'waveform_metrics','PCA_features','acg_metrics','monoSynaptic_connections','theta_metrics','spatial_metrics','event_metrics','psth_metrics'},@iscell);

addParameter(p,'keepCellClassification',true,@islogical);
addParameter(p,'manualAdjustMonoSyn',true,@islogical);
addParameter(p,'getWaveformsFromDat',true,@islogical);
addParameter(p,'includeInhibitoryConnections',false,@islogical);
addParameter(p,'showGUI',false,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'forceReloadSpikes',false,@islogical);
addParameter(p,'forceReloadWaveformBasedMetrics',false,@islogical);
addParameter(p,'submitToDatabase',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveAs','cell_metrics',@isstr);
addParameter(p,'saveBackup',true,@islogical);
addParameter(p,'fileFormat','mat',@isstr);
addParameter(p,'transferFilesFromClusterpath',true,@islogical);

% Plot related parameters
addParameter(p,'showFigures',false,@islogical);
addParameter(p,'showWaveforms',true,@islogical);
addParameter(p,'summaryFigures',false,@islogical);
addParameter(p,'debugMode',false,@islogical);     
 
parse(p,varargin{:})

sessionID = p.Results.sessionID;
sessionin = p.Results.sessionName;
sessionStruct = p.Results.session;
basepath = p.Results.basepath;
parameters = p.Results;
timerCalcMetrics = tic;

% Verifying required toolboxes are installed
installedToolboxes = ver;
installedToolboxes = {installedToolboxes.Name};
requiredToolboxes = {'Curve Fitting Toolbox','Signal Processing Toolbox','Statistics and Machine Learning Toolbox'};

missingToolboxes = requiredToolboxes(~ismember(requiredToolboxes,installedToolboxes));
if ~isempty(missingToolboxes)
    for i = 1:numel(missingToolboxes)
        warning(['A toolbox required by CellExplorer must be installed: ' missingToolboxes{i}]);
    end
end

% Load MergePoints
targetFile = dir('*MergePoints.events.mat');
if ~isempty(targetFile)
    load(targetFile.name);
else
    error('MergePoints not found. Quitting...');
end

excludeManipulationIntervals = [];
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
    
for ii = 1:length(MergePoints.foldernames)
    ts = MergePoints.timestamps(ii,:);
    
    cell_metrics = ProcessCellMetrics('restrictToIntervals',ts,'session',parameters.session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',parameters.excludeMetrics,'forceReload',parameters.forceReload,'saveMat',false);
    
end


end