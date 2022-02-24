function  preprocessSession(varargin)

%         bz_PreprocessSession(varargin)

%   Master function to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones.
% 

% INPUTS
%   <options>         optional list of property-value pairs (see table below)
%   basepath          Basepath for experiment. It contains all session
%                       folders. If not provided takes pwd.
%   analogCh          List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
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
%   sessionSummary    Default, 'all'. 
%  tracking_pixel_cm  Default, 0.1149
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
addParameter(p,'basepath',pwd,@isdir); % by default, current folder
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'cleanArtifacts',false);
addParameter(p,'medianSubstr',true);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'sessionSummary',true,@islogical);

% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
analogCh = p.Results.analogCh;
spikeSort = p.Results.spikeSort;
getPos = p.Results.getPos;
cleanArtifacts = p.Results.cleanArtifacts;
medianSubstr = p.Results.medianSubstr;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
sessionSummary = p.Results.sessionSummary;

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
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),[basename '.xml']);
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));
    end
    cd(basepath);
end

%% Loading metadata
try
    session = sessionTemplate(pwd,'showGUI',false); % 
    save([basename '.session.mat'],'session');
catch
    warning('it seems that CellExplorer is not on your path');
end

%% Concatenate sessions
cd(basepath);
disp('Concatenate session folders...'); 
concatenateDats(pwd,0,1);

%% Get analog and digital pulses
if  ~isempty(analogCh)
    [pulses] = getAnalogPulses('analogCh',analogCh,'manualThr', false);
end
if ~isempty(dir('*digitalIn.mat'))
    digitalIn = getDigitalIn('all','fs',session.extracellular.sr); 
end

%% Remove stimulation artifacts
if iscell(cleanArtifacts) || cleanArtifacts
    if iscell(cleanArtifacts)
        pulArtifacts_analog = getAnalogPulses('analogCh',cleanArtifacts{1});
        pulArtifacts_dig = [];
        for ii = cleanArtifacts{2}
            disp(ii);
            pulArtifacts_dig = [pulArtifacts_dig; digitalIn.ints{ii}(:)];
        end
        pulArtifacts = sort([pulArtifacts_analog.timestamps(:); pulArtifacts_dig]);
    else
        pulArtifacts = pulses.timestamps(:);
    end
    cleanPulses(pulArtifacts);
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
cd(basepath);

%% Get tracking positions 
if getPos
    getSessionTracking('convFact',tracking_pixel_cm,'roiTracking','manual'); 
end

if sessionSummary
    computeSessionSummary;
end

end

