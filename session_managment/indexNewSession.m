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
addParameter(p,'rippleChannel',[],@isnumeric);% manually selecting ripple Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'SWChannel',[],@isnumeric); % manually selecting SW Channel in case getHippocampalLayers does not provide a right output
addParameter(p,'digitalChannelsList',[],@isnumeric);
addParameter(p,'analogChannelsList',[],@isnumeric);
addParameter(p,'promt_hippo_layers',false,@islogical);
addParameter(p,'manual_analog_pulses_threshold',false,@islogical);
addParameter(p,'removeDatFiles',true,@islogical);
addParameter(p,'removeDat',false,@islogical);
addParameter(p,'copyFiles',true,@islogical);
addParameter(p,'copyPath',[],@isdir);
addParameter(p,'bazler_ttl_channel',[],@isnumeric);
addParameter(p,'forceAnalogPulses',false,@islogical);
addParameter(p,'forceDigitalPulses',false,@islogical);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'probe_type',[],@ischar); % 

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
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
digitalChannelsList = p.Results.digitalChannelsList;
analogChannelsList = p.Results.analogChannelsList;
promt_hippo_layers = p.Results.promt_hippo_layers;
manual_analog_pulses_threshold = p.Results.manual_analog_pulses_threshold;
removeDatFiles = p.Results.removeDatFiles;
removeDat = p.Results.removeDat;
copyFiles = p.Results.copyFiles;
copyPath = p.Results.copyPath;
bazler_ttl_channel = p.Results.bazler_ttl_channel;
forceAnalogPulses = p.Results.forceAnalogPulses;
forceDigitalPulses = p.Results.forceDigitalPulses;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
probe_type = p.Results.probe_type;

%% Creates a pointer to the folder where the index variable is located
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.csv']));
    if isempty(indexedProjects_path)
        disp('No indexed Projects .csv file found. Lets create one !' );
        directory = what(hippoCookBook_path);
        cd(directory.path);
        allSessions = [];
        save([indexedProjects_name,'.csv'],'allSessions');
        indexedProjects_path = fileparts(which([indexedProjects_name,'.csv']));
    end
end

cd(basepath)

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
    
    if ~isfield(session, 'analysisTags')
        session.analysisTags = [];
    end
    if ~isfield(session.analysisTags,'digital_optogenetic_channels')
        session.analysisTags.digital_optogenetic_channels = digitalChannelsList;
    end
    if ~isfield(session.analysisTags,'analog_optogenetic_channels')
        session.analysisTags.analog_optogenetic_channels = analogChannelsList;
    end
    if isempty(rejectChannels)
        rejectChannels = session.channelTags.Bad.channels; % 1-index
    end
    if ~isfield(session.analysisTags,'bazler_ttl_channel')
        session.analysisTags.bazler_ttl_channel = bazler_ttl_channel;
    end
    save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
catch
    warning('it seems that CellExplorer is not on your path');
end

if ~isempty(probe_type)
    switch lower(probe_type)
        case  {'fuckyou',lower('A5x12-16-Buz-lin-5mm-100-200-160-177')}
            disp('Probe founded!!');
            directory = what(hippoCookBook_path);
            coord_path = dir([directory.path filesep 'session_files' filesep 'probes_coordinates' filesep 'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.chanCoords.channelInfo.mat']);
            load([coord_path.folder filesep coord_path.name],'chanCoords');
            save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords');
            
            figure
            hold on
            plot(chanCoords.x, chanCoords.y,'.','MarkerSize',10,'color',[.8 .8 .8]);
            for ii = 1:length(chanCoords.x)
                text(chanCoords.x(ii)+1, chanCoords.y(ii)-1, num2str(ii));
            end
            xlabel('um'); ylabel('um');
            title(['A5x12-16-Buz-lin-5mm-100-200-160-177 probe'],'FontWeight','normal');
            xlim([min(chanCoords.x)-100 max(chanCoords.x)+100]);
            ylim([min(chanCoords.y)-100 max(chanCoords.y)+100]);
            set(gca,'TickDir','out');
            mkdir('SummaryFigures');
            saveas(gcf,['SummaryFigures\probe_layout.png']);
            
            session = loadSession(basepath);
            session.animal.probeImplants{1}.probe = 'A5x12-16-Buz-lin-5mm-100-200-160-177';
            session.animal.probeImplants{1}.layout = 'A5x12-16-Buz-lin-5mm-100-200-160-177';
            session.extracellular.chanCoords = chanCoords;
            save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
        otherwise
            disp('Probe not supported yet...');
    end
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
if forceAnalogPulses || isempty(dir([session.general.name,'_original.dat']))
    disp('Getting analog Pulses...')
    pulses = getAnalogPulses('analogChannelsList',analogChannelsList,'manualThr',manual_analog_pulses_threshold,'overwrite',force_analogPulsesDetection); % 1-index
else
    try
        if ~isempty(dir([session.general.name,'.pulses.events.mat']))
            file = dir([session.general.name,'.pulses.events.mat']);
            load(file.name);
        end
    catch
        warning('Problems with analogPulses');
    end
end
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
getAverageCCG('force',true);

% pulses.analogChannel = analogCh;
% save([session.general.name,'.pulses.events.mat'],'pulses');
optogeneticResponses = getOptogeneticResponse('numRep',500,'force',true,'duration_round_decimal',1);

%% 8. Check Brain Events
% Trying changes in detecUD_temp
% 8.1 Up and downs
UDStates = detectUD('plotOpt', true,'forceDetect',true','NREMInts','all');
psthUD = spikesPsth([],'eventType','slowOscillations','numRep',500,'force',true);

% 8.2 Ripples
ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true);
psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true);

% 8.3 Theta intervals
thetaEpochs = detectThetaEpochs();

%% 9. Phase Modulation
% LFP-spikes modulation
[rippleMod,SWMod,thetaMod,lgammaMod,hgammaMod] = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel);

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
    getSessionTracking('convFact',tracking_pixel_cm,'roiTracking','manual');
    getSessionArmChoice('task','alternation');
    behaviour = getSessionLinearize('forceReload',false);  
    firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
    placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %% ,'maxSize',.75,'sepEdge',0.03
    firingTrialsMap = firingMapPerTrial('force',true);
    spatialModulation = getSpatialModulation('force',true);
catch
    warning('Not possible to run spatial modulation...');
end

%% 12. Events modulation
try 
    behaviour = getSessionLinearize;
    psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',100,'saveMat',false,...
        'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
    psth_rReward = spikesPsth([behaviour.events.rReward],'numRep',100,'saveMat',false,...
        'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
    psth_reward = spikesPsth([behaviour.events.lReward; behaviour.events.rReward],'numRep',100,'saveMat',false,...
        'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
    
    psth_intersection = spikesPsth([behaviour.events.intersection],'numRep',100,'saveMat',false,...
        'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);
    psth_startPoint = spikesPsth([behaviour.events.startPoint],'numRep',100,'saveMat',false,...
        'min_pulsesNumber',5,'winSize',6,'event_ints',[0 0.2],'winSizePlot',[-2 2],'binSize',0.01, 'win_Z',[-3 -1]);

    behaviour.psth_lReward = psth_lReward;
    behaviour.psth_rReward = psth_rReward;
    behaviour.psth_reward = psth_reward;
    behaviour.psth_intersection = psth_intersection;
    behaviour.psth_startPoint = psth_startPoint; 
    behavior = behaviour; % british to american :)
    save([basenameFromBasepath(pwd) '.behavior.cellinfo.mat'],'behavior');
end

%% 13. Speed Score

try
    speedCorr = getSpeedCorr('numQuantiles',20);
end

%% 14. Summary per cell
getACGPeak;
getSummaryPerCell;

%% 14. Indexing
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
generalPath = [session.animal.name,'\',session.general.name];
sessionName = session.general.name;

project = session.general.projects;
% updated indexedSession table
sessionsTable = readtable([indexedProjects_path filesep indexedProjects_name,'.csv']); % the variable is called allSessions
% new table entry

optogenetics = cell(0);
for ii = 1:length(session.animal.opticFiberImplants)
    optogenetics{1, length(optogenetics)+1} = session.animal.opticFiberImplants{ii}.opticFiber;
    optogenetics{1, length(optogenetics)+1} = ' ';
end
optogenetics(end) = [];

behav = cell(0); 
for i = 1:length(session.epochs)
    if strcmpi(session.epochs{i}.behavioralParadigm, 'Maze')
        behav{1, length(behav)+1} = lower(session.epochs{i}.environment);
        behav{1, length(behav)+1} = ' ';
    end
end

if ~isempty(behav)
    behav(end) = [];
    if isempty(behav)
        behav{1,1} = 'no';
    end
else
    behav{1,1} = 'no';
end

spikes = loadSpikes;

fn = fieldnames(session.brainRegions);
brainRegions = cell(0);
for jj = 1:length(fn)
    brainRegions{1,length(brainRegions)+1} = fn{jj};
    brainRegions{1,length(brainRegions)+1} = ' ';
end    
brainRegions(end) = [];

sessionEntry = {lower(sessionName), lower(session.animal.name), lower(generalPath), lower(session.animal.strain),...
    lower(session.animal.geneticLine), lower([optogenetics{:}]), [behav{:}], spikes.numcells,  [brainRegions{:}], project};
sessionEntry = cell2table(sessionEntry,"VariableNames",["SessionName", "Subject", "Path", "Strain", "GeneticLine", "Optogenetics", "Behavior", "numCells", "brainRegions", "Project"]);
sessionsTable = [sessionsTable; sessionEntry];
writetable(sessionsTable,[indexedProjects_path filesep indexedProjects_name,'.csv']); % the variable is called allSessions


% Lets do a push for git repository
cd(indexedProjects_path);
% Git add variable to the repository
commandToExecute = ['git add ', indexedProjects_name,'.csv']
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
            try
                cd(MergePoints.foldernames{i})
                if ~isempty(dir('amplifier*.dat'))
                    file = dir('amplifier*.dat');
                    delete(file.name);
                end
                cd(basepath)
            catch
                warning('Problem removing session folders...');
            end
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
