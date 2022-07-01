function  preprocessSession_pablo(varargin)

%         bz_PreprocessSession(varargin)

%   Master function to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones.
% 

% INPUTS
%   <options>         optional list of property-value pairs (see table below)
%   basepath          Basepath for experiment. It contains all session
%                       folders. If not provided takes pwd.
%   analysisPath        local path to run preprocessing. If empty,
%                       preprocessing will be done in basepath, if any
%                       folder selected, files will be copied from basepath
%                       to analysisPath, run the analysis to speed up and
%                       then copy back again to basepath and delete
%                       analysisPath
%   analogChannelsList  
%                     List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum          Force make folder summary (overwrite, if necessary). Default false.
%   cleanArtifacts    Remove artifacts from dat file (false by default). If
%                       true, remove artifacts from all Analog events. It also
%                       accepts a two rows cell with the analog channel
%                       (cleanArtifacts{1}) and the digital channels
%                       (cleanArtifacts{2}) to be used.
%   stateScore        Run automatic brain state detection with SleepScoreMaster. Default false.
%   spikeSort         Run automatic spike sorting using Kilosort. Default true.
%   getPos            Get tracking positions. Default true. 
%   medianSubstr      Perform median substraction in dat file before
%                       kilosort. Careful!! it would compromises dat file!
%                       (default false). If scalar, perform median
%                       substraction in those channels.
%   sessionSummary    Default, 'false'. 
%   tracking_pixel_cm Default, 0.1149
%   digitalChannelsList     Array of channel to perform 'digitalPulses'
%                     if summary is done, otherwise [] 
%   bazler_ttl_channel Channel for bazler tracking ttl
%
%  HISTORY: 
%     - Created based on sessionsPipeline
%     - Updated for hippoCookBook, 2022
%
%  TO DO:
%   - Include Kilosort2 support
%   - Improve auto-clustering routine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'analysisPath',[]); % Local paht to run the anaysis, is em
addParameter(p,'analogChannelsList',[],@isnumeric);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'cleanArtifacts',false);
addParameter(p,'medianSubstr',true);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'sessionSummary',true,@islogical);
addParameter(p,'digitalChannelsList',[],@isnumeric);
addParameter(p,'manualThr',true,@islogical);
addParameter(p,'bazler_ttl_channel',10,@isnumeric);
addParameter(p,'anymaze_ttl_channel',[],@isnumeric);
addParameter(p,'getDigitalInputBySubfolders',true,@islogical);
addParameter(p,'anyMaze',true,@islogical);

% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
analysisPath = p.Results.analysisPath;
analogChannelsList = p.Results.analogChannelsList;
spikeSort = p.Results.spikeSort;
getPos = p.Results.getPos;
cleanArtifacts = p.Results.cleanArtifacts;
medianSubstr = p.Results.medianSubstr;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
sessionSummary = p.Results.sessionSummary;
digitalChannelsList = p.Results.digitalChannelsList;
manualThr = p.Results.manualThr;
bazler_ttl_channel = p.Results.bazler_ttl_channel;
anymaze_ttl_channel = p.Results.anymaze_ttl_channel;
getDigitalInputBySubfolders = p.Results.getDigitalInputBySubfolders;
anyMaze = p.Results.anyMaze;

if ~exist('basepath') || isempty(basepath)
    basepath = uigetdir; % select folder
end
cd(basepath);

%% deals with xml.
if strcmp(basepath(end),filesep)
    basepath = basepath(1:end-1);
end
[~,basename] = fileparts(basepath);

disp('Check xml...');
if isempty(dir([basename '.xml'])) && isempty(dir('global.xml'))
    disp('No xml global file! Looking for it...');
    allpath = strsplit(genpath(basepath),';'); % all folders
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*amplifier*.xml');
        ii = ii + 1;
    end    
    if isempty(xmlFile)
        % Looking global.xml in general folder
        cd(basepath)
        cd ..
        xmlFile = dir('*global*.xml');
    end 
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),[basepath,filesep,basename '.xml']);
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));
    end
    cd(basepath);
end

%% Loading metadata
try
    session = sessionTemplate(pwd,'showGUI',false); % 
    session.channels = 1:session.extracellular.nChannels;
    session.analysisTags.digital_optogenetic_channels = digitalChannelsList;
    session.analysisTags.analog_optogenetic_channels = analogChannelsList;
    if ~isempty(bazler_ttl_channel)
        session.analysisTags.bazler_ttl_channel = bazler_ttl_channel;
    end
    if ~isempty(anymaze_ttl_channel)
        session.analysisTags.anymaze_ttl_channel = anymaze_ttl_channel;
    end
    save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
catch
    warning('it seems that CellExplorer is not on your path');
end

if ~isempty(analysisPath)
    cd(analysisPath);
    mkdir(session.general.name);
    copyfile(basepath,[analysisPath,'\',session.general.name]);
    disp('Copying files to analysis path folder...');
    cd([analysisPath,'\',session.general.name]);
    disp('Copied files. Peforming preprocessSession...')
end
    
%% Concatenate sessions
disp('Concatenate session folders...'); 
concatenateDats(pwd,0,1);

%% Get analog and digital pulses
if  ~isempty(analogChannelsList)
    try
        [pulses] = getAnalogPulses('analogChannelsList',analogChannelsList,'manualThr', manualThr);
    catch
        warning('No analog pulses detected');
    end
end
if ~isempty(dir('*digitalIn.dat')) 
    digitalIn = getDigitalIn('all','fs',session.extracellular.sr); 
end

if getDigitalInputBySubfolders
    try
        digitalIn = getDigitalInBySubfolders('all','fs',session.extracellular.sr);
    end
end

% digitalIn = pap_getDigitalIn('all','fs',session.extracellular.sr);

%% Remove stimulation artifacts
try
    if iscell(cleanArtifacts) || cleanArtifacts
        if iscell(cleanArtifacts)
            pulArtifacts_analog = getAnalogPulses('analogChannelsList',cleanArtifacts{1});
            if isempty(pulArtifacts_analog)
                pulArtifacts_analog.timestamps = [];
            end            
            
            pulArtifacts_dig = [];
          digitalIn = getDigitalIn;
            if isempty(digitalIn) || isempty(cleanArtifacts{2})
            else
                for ii = cleanArtifacts{2}
                    pulArtifacts_dig = [pulArtifacts_dig; digitalIn.ints{ii}];
                end
            end
            pulArtifacts = [pulArtifacts_analog.timestamps; pulArtifacts_dig];
        else
            pulArtifacts = pulses.timestamps;
        end
        removeStimulationArtifacts(pulArtifacts);
    end
catch
    warning('remove stimulation artifacts not possible. Probably digital input is empty.')
end

%% Make LFP
if isempty(dir('*.lfp'))
    disp('Creating .lfp file. This could take a while...');
    ResampleBinary(strcat(basename,'.dat'),strcat(basename,'.lfp'),...
        session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
end

%% MEDIAN SUBS
if isempty(dir([session.general.name '_original.dat']))
    if islogical(medianSubstr) && medianSubstr
        medianSubtraction(pwd);
    elseif medianSubstr
        medianSubtraction(pwd,'ch',medianSubstr);
    end
else
    warning('Session was already median-subtracted. Spiking...');
end

%% Kilosort concatenated sessions
if spikeSort
    if  isempty(dir('*Kilosort*')) % if not kilosorted yet
        fprintf(' ** Kilosorting session...');
        
        KiloSortWrapper;
        kilosortFolder = dir('*Kilosort*');
        try PhyAutoClustering(strcat(kilosortFolder.folder,filesep,kilosortFolder.name)); % autoclustering
        catch
            warning('PhyAutoClustering not possible!!');
        end
        if exist('phyLink') && ~isempty(phyLink) % move phy link to
            kilosort_path = dir('*Kilosort*');
            try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
            end
        end        
    else 
        warning('Session was already run kilosort in. Spiking...');
    end
end

if ~isempty(analysisPath)
    cd([analysisPath,'\',session.general.name])
else
    cd(basepath)
end

%% Get tracking positions 
if getPos
    getSessionTracking('convFact',tracking_pixel_cm,'roiTracking','manual','anyMaze',anyMaze); 
end

if sessionSummary
    cd(basepath);
    session = sessionTemplate(pwd,'showGUI',false);
    save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
    computeSessionSummary_pablo('digitalChannelsList',digitalChannelsList,'analogChannelsList',analogChannelsList);
end

if ~isempty(analysisPath)
    [success] = copyfile([analysisPath,'\',session.general.name],basepath);
    if success
        cd(analysisPath);
        try
            rmdir([analysisPath,'\',session.general.name],'s')
        catch
            disp('Not possible to remove folder...');
            fclose('all');
            CloseNoPrompt;
            rmdir([analysisPath,'\',session.general.name],'s');
            disp('Folder deleted');
        end
    end
end

function CloseNoPrompt
%Close all editor windows without prompting
%Active Editor;
hEditor = matlab.desktop.editor.getActive;
%Close all files.
while ~isempty(hEditor)
  closeNoPrompt(hEditor);
  hEditor = matlab.desktop.editor.getActive;
end
