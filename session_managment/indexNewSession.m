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

%% Pablo Abad and Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'rejectChannels',[],@isnumeric); % 0-index
addParameter(p,'project','Undefined',@isstring);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'force_analogPulsesDetection',true,@islogical);
addParameter(p,'force_loadingSpikes',true,@islogical);
addParameter(p,'excludeManipulationIntervals',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric); % manually selecting SW Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'digitalChannelsList',[],@isnumeric);
addParameter(p,'analogChannelsList',[],@isnumeric);
addParameter(p,'promt_hippo_layers',false,@islogical);
addParameter(p,'manual_analog_pulses_threshold',false,@islogical);
addParameter(p,'removeDatFiles',true,@islogical);
addParameter(p,'removeDat',false,@islogical);
addParameter(p,'copyFiles',true,@islogical);
addParameter(p,'copyPath',[],@isdir);

parse(p,varargin{:})

basepath = p.Results.basepath;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;
rejectChannels = p.Results.rejectChannels;
project = p.Results.project;
indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
hippoCookBook_path = p.Results.hippoCookBook_path;
force_analogPulsesDetection = p.Results.force_analogPulsesDetection;
force_loadingSpikes = p.Results.force_loadingSpikes;
excludeManipulationIntervals = p.Results.excludeManipulationIntervals;
SWChannel = p.Results.SWChannel;
digitalChannelsList = p.Results.digitalChannelsList;
analogChannelsList = p.Results.analogChannelsList;
promt_hippo_layers = p.Results.promt_hippo_layers;
manual_analog_pulses_threshold = p.Results.manual_analog_pulses_threshold;
removeDatFiles = p.Results.removeDatFiles;
removeDat = p.Results.removeDat;
copyFiles = p.Results.copyFiles;
copyPath = p.Results.copyPath;

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
keyboard

%% By default looks for Synology and copy files to it, it specified copy files to the specified folder
if isempty(copyPath)
    % Let's find the packrat synology 
    F = getdrives();
    for i = 1:length(F)
        if strcmpi(driveName(cell2mat(strsplit(F{i},[':',filesep]))),'packrat')
            copyPath = [F{i},'data'];
            cd(copyPath);
        end
    end
end
cd(basepath)
%% 1. Runs sessionTemplate
try
    session = loadSession(basepath);
    session.channels = 1:session.extracellular.nChannels;
    if ~isfield(session.analysisTags,'digital_optogenetic_channels')
        session.analysisTags.digital_optogenetic_channels = digitalChannelsList;
    end
    if ~isfield(session.analysisTags,'analog_optogenetic_channels')
        session.analysisTags.analog_optogenetic_channels = analogChannelsList;
    end
    if isempty(rejectChannels)
        rejectChannels = session.channelTags.Bad.channels; % 1-index
    end
    save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
catch
    warning('it seems that CellExplorer is not on your path');
end
session = sessionTemplate(basepath,'showGUI',true);
%% 2. Remove previous cellinfo.spikes.mat and computes spikes again (manual clustered)
disp('Loading Spikes...')
if force_loadingSpikes
    if ~isempty(dir([session.general.name ,'.spikes.cellinfo.mat']))
        file = dir([session.general.name ,'.spikes.cellinfo.mat']);
        delete(file.name);
    end
end   
spikes = loadSpikes('forceReload',force_loadingSpikes);

%% 3. Analog pulses detection
disp('Getting analog Pulses...')
pulses = getAnalogPulses('analogChannelsList',analogChannelsList,'manualThr',manual_analog_pulses_threshold,'overwrite',force_analogPulsesDetection); % 1-index

%% 4. Check Sleep Score
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods, 'overwrite', true);
bz_ThetaStates(pwd);

%% 5. Power Profiles
powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'forceDetect',true);
powerProfile_gamma = powerSpectrumProfile(gamma_bandpass,'showfig',true,'forceDetect',true);
powerProfile_hfo = powerSpectrumProfile(hfo_bandpass,'showfig',true,'forceDetect',true);

%% 6. Getting Hippocampal Layers
[hippocampalLayers] = getHippocampalLayers('force',true,'promt',promt_hippo_layers);

%% 7. Spike Features
spikeFeatures;
getAverageCCG;
% pulses.analogChannel = analogCh;
% save([session.general.name,'.pulses.events.mat'],'pulses');
optogeneticResponses = getOptogeneticResponse('numRep',500,'force',true);

%% 8. Check Brain Events
% Trying changes in detecUD_temp
% 8.1 Up and downs
UDStates = detectUD('plotOpt', true,'forceDetect',true','NREMInts','all');
psthUD = spikesPsth([],'eventType','slowOscillations','numRep',500);

% 8.2 Ripples
ripples = rippleMasterDetector('SWChannel',SWChannel,'force',true);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500);

% 8.3 Theta intervals
thetaEpochs = detectThetaEpochs;

%% 9. Phase Modulation
% LFP-spikes modulation
[rippleMod,SWMod,thetaMod,lgammaMod,hgammaMod] = computePhaseModulation('SWChannel',SWChannel);

%% 10. Cell metrics
% Exclude manipulation intervals for computing CellMetrics
try
    if ~isempty(dir([session.general.name,'.optogeneticPulses.events.mat']))
        file = dir([session.general.name,'.optogeneticPulses.events.mat']);
        load(file.name);
    end
        excludeManipulationIntervals = optoPulses.stimulationEpochs;
catch
    warning('Not possible to get manipulation periods. Running CellMetrics withouth excluding manipulation epochs');
end
cell_metrics = ProcessCellMetrics('session', session,'excludeIntervals',excludeManipulationIntervals,'excludeMetrics',{'deepSuperficial'});

%% 11. Spatial modulation
try
    behaviour = getSessionLinearize('forceReload',false);  
    firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',false);
    placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
    firingTrialsMap = firingMapPerTrial;
catch
    warning('Not possible to run spatial modulation...');
end

%% 11. Indexing
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
generalPath = [session.animal.name,'\',session.general.name];
sessionName = session.general.name;
load([indexedProjects_path filesep indexedProjects_name,'.mat']); % the variable is called allSessions
allSessions.(sessionName).path = generalPath;
allSessions.(sessionName).name = session.animal.name;
allSessions.(sessionName).strain = session.animal.strain;
allSessions.(sessionName).geneticLine = session.animal.geneticLine;
allSessions.(sessionName).brainRegions = session.brainRegions;
if isfield(session.animal,'opticFiberImplants')
    for i = 1:length(session.animal.opticFiberImplants)
        allSessions.(sessionName).optogenetics{i} = session.animal.opticFiberImplants{i}.opticFiber;
    end
else
    allSessions.(sessionName).optogenetics = NaN;
end
behav = NaN;
for i = 1:length(session.epochs)
    if strcmpi(session.epochs{i}.behavioralParadigm, 'Maze')
        behav = [behav session.epochs{i}.environment];
    end
end
allSessions.(sessionName).behav = behav;
allSessions.(sessionName).project = session.general.projects;
save([indexedProjects_path filesep indexedProjects_name,'.mat'],'allSessions');

% Lets do a push for git repository
cd(indexedProjects_path);
% Git add variable to the repository
commandToExecute = ['git add ', indexedProjects_name,'.mat']
system(commandToExecute);
% Git Commit
commentToCommit = ['Added Session: ' session.general.name];
commandToExecute = ['git commit -m "' commentToCommit '"'];
system(commandToExecute);
% Git Push
commandToExecute = ['git push'];
system(commandToExecute);

cd(basepath)

%% Removing dat files before copying files to buzsakilab or synology
if removeDatFiles
    % Remove _original and _temp .dat
    if ~isempty(dir([session.general.name,'_original.dat']))
        delete([session.general.name,'_original.dat']);
    end
    if ~isempty(dir([session.general.name,'_temp.dat']))
        delete([session.general.name,'_temp.dat']);
    end
    
    % Remove amplifier*.dat in subfolders
    if ~isempty(dir([session.general.name,'.MergePoints.events.mat']))
        file = dir([session.general.name,'.MergePoints.events.mat']);
        load(file.name)
        
        for i = 1:length(MergePoints.foldernames)
            cd(MergePoints.foldernames{i})
            if ~isempty(dir('amplifier*.dat'))
                file = dir('amplifier*.dat');
                delete(file.name);
            end
            cd(basepath)
        end
    end
    
    % Remove kilosort .phy
    if ~isempty(dir('Kilosort*'))
        file = dir('Kilosort*');
        cd(file.name);
        if exist('.phy','dir')
            rmdir('.phy','s');
        end
        cd(basepath);
    end
end

if removeDat
    if ~isempty(dir([session.general.name,'.dat']))
        file = dir([session.general.name,'.dat']);
        delete(file.name);
    end
end

%% TO DO. Copy files to remote and delete the folder in this computer
if copyFiles
    if ~exist([copyPath,'\',session.animal.name],'dir')
        mkdir([copyPath,'\',session.animal.name]);
    end
    [success] = copyfile(session.general.basePath,[copyPath,'\',session.animal.name,'\',session.general.name]);
    if success
        cd ..
        rmdir(session.general.basePath,'s');
    end
end

end
