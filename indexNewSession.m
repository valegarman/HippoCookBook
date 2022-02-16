function [] = indexNewSession(varargin)

%       [] = indexSession_InterneuronsLibrary(varargin)

% This function runs standard preAnalysis
% Based on index_session_script_InterneuronsLibrary by MV 2020
% 1. Runs sessionTemplate
% 2. remove previous cellinfo.spikes.mat and computes spikes again (
%       manual clustered)
% 3. remove previous opto pulses file (pulses.events.mat) and re-runs opto
%       pulses analysis
% 4. Runs and Check SleepScore
% 5. Power Profiles
% 6. check UD events
% 7. Ripples analysis
% 8. Cell Metrics and CellExplorer
% 9. Spikes features
% 10. Saving index path

% TO DO:
% Deactivate by inputs loadSpikes or other analysis

%% Pablo Abad 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'analogCh',65,@isnumeric); % 0-index
addParameter(p,'rejectChannels',[],@isnumeric); % 0-index
addParameter(p,'project','Undefined',@isstring);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'force_analogPulsesDetection',true,@islogical);
addParameter(p,'excludeManipulationIntervals',[],@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;
analogCh = p.Results.analogCh;
rejectChannels = p.Results.rejectChannels;
project = p.Results.project;
indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
hippoCookBook_path = p.Results.hippoCookBook_path;
force_analogPulsesDetection = p.Results.force_analogPulsesDetection;
excludeManipulationIntervals = p.Results.excludeManipulationIntervals;

keyboard;
%% Creates a pointer to the folder where the index variable is located
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .mat variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.mat']));
    if isempty(indexedProjects_path)
        disp('No indexed Projects .mat file found. Lets create one !' );
        directory = what(hippoCookBook_path);
        cd(directory.path);
        allSessions = [];
        save([indexedProjects_name,'.mat'],'allSessions');
        indexedProjects_path = fileparts(which([indexedProjects_name,'.mat']));
    end
end

cd(basepath)

%% 1. Runs sessionTemplate
session = sessionTemplate(basepath,'showGUI',true);
% Parsing rejectChannels from session.mat file in case rejectChannels is empty
if isempty(rejectChannels)
    rejectChannels = session.channelTags.Bad.channels; % 1-index
end
% Creating a field in session.mat called channels (1-index)
session.channels = 1:session.extracellular.nChannels;
save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');

%% 2. Remove previous cellinfo.spikes.mat and computes spikes again (manual clustered)

if ~isempty(dir([basepath filesep session.general.name ,'.spikes.cellinfo.mat']))
    disp('Loading and deleting old spikes.cellinfo.mat file ...');
    file = dir([basepath filesep session.general.name ,'.spikes.cellinfo.mat']);
    delete(file.name);
else
    ('spikes.cellinfo.mat does not exist !');
end

disp('Loading Spikes...')
spikes = loadSpikes;


%% 3. Analog pulses detection

disp('Getting analog Pulses...')
pulses = getAnalogPulses('analogCh',analogCh,'manualThr',false,'overwrite',force_analogPulsesDetection); % 1-index

%% 4. Check Sleep Score
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods, 'overwrite', true);
bz_ThetaStates(pwd);

%% 5. Power Profiles
powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'forceDetect',true);
powerProfile_gamma = powerSpectrumProfile(gamma_bandpass,'showfig',true,'forceDetect',true);
powerProfile_hfo = powerSpectrumProfile(hfo_bandpass,'showfig',true,'forceDetect',true);

%% 6. Getting Hippocampal Layers
[hippocampalLayers] = getHippocampalLayers('force',true);

%% 7. Check Brain Events
% Trying changes in detecUD_temp
% 7.1 Up and downs
UDStates = detectUD('plotOpt', true,'forceDetect',true','NREMInts','all');

% 7.2 Ripples
ripples = rippleMasterDetector('SWChannel',48);

% 7.3 Theta intervals
thetaEpochs = detectThetaEpochs;

%% 10. Cell metrics
% Exclude manipulation intervals for computing CellMetrics
try
    excludeManipulationIntervals = pulses.intsPeriods;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});
%% 8. Cell metrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});
% cell_metrics = CellExplorer('metrics',cell_metrics);

%% 9. Spike Features
spikeFeatures()
% pulses.analogChannel = analogCh;
% save([session.general.name,'.pulses.events.mat'],'pulses');
optogeneticResponses = getOptogeneticResponse('numRep',100);

%% 10. Indexing
session = sessionTemplate(basepath,'showGUI',false);
currentPath = split(pwd,':'); currentPath = currentPath{end};
sessionName = session.general.name;
load([indexedProjects_path filesep indexedProjects_name,'.mat']); % the variable is called allSessions
allSessions.(sessionName).path = currentPath;
allSessions.(sessionName).name = session.animal.name;
allSessions.(sessionName).strain = session.animal.strain;
allSessions.(sessionName).geneticLine = session.animal.geneticLine;
if isfield(session.animal,'opticFiberImplants')
    for i = 1:length(session.animal.opticFiberImplants)
        allSessions.(sessionName).optogenetics{i} = session.animal.opticFiberImplants{i}.opticFiber;
    end
end
behav = [];
for i = 1:length(session.epochs)
    if strcmpi(session.epochs{i}.behavioralParadigm, 'Maze')
        behav = [behav session.epochs{i}.environment];
    end
end
allSessions.(sessionName).behav = behav;
allSessions.(sessionName).project = sessions.general.projects;
allSessions.(sessionName).tag = 1;
save([indexedProjects_path filesep indexedProjects_name,'.mat'],'allSessions');

% Lets do a push for git repository
cd(indexedProjects_path)
% Git add variable to the repository
commandToExecute = ['git add ', indexedProjects_name,'.mat']
system(commandToExecute)
% Git Commit
commentToCommit = ['Added Session: ' session.general.name];
commandToExecute = ['git commit -m "' commentToCommit '"'];
system(commandToExecute)
% Git Push
commandToExecute = ['git push'];
system(commandToExecute)

end
