function  batch_preprocessSession(varargin)

%   Runs the master function preprocessSession throughout all sessions in
%   the basepath without a kilosort folder on them.
% 

% INPUTS
%   <options>         optional list of property-value pairs (see table below)
%   basepath          Basepath for experiment. It contains all session
%                       folders. If not provided takes pwd.
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
%   tracking_pixel_cm  Default, 0.1149
%   digitalChannelsList     Array of channel to perform 'digitalPulses'
%                     if summary is done, otherwise [] 
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
addParameter(p,'analogChannelsList',[],@isnumeric);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',false,@islogical);
addParameter(p,'cleanArtifacts',true);
addParameter(p,'medianSubstr',true);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'sessionSummary',false,@islogical);
addParameter(p,'digitalChannelsList',[],@isnumeric);

% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
analogChannelsList = p.Results.analogChannelsList;
spikeSort = p.Results.spikeSort;
getPos = p.Results.getPos;
cleanArtifacts = p.Results.cleanArtifacts;
medianSubstr = p.Results.medianSubstr;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
sessionSummary = p.Results.sessionSummary;
digitalChannelsList = p.Results.digitalChannelsList;

% batch processing...
all_folders = dir(basepath);
for ii = 1:size(all_folders,1)
    if all_folders(ii).isdir && ~strcmpi(all_folders(ii).name,'.') && ~strcmpi(all_folders(ii).name,'..')
        cd([all_folders(ii).folder filesep all_folders(ii).name]);
        kilosortFolder = dir('*Kilosort*');
        if isempty(kilosortFolder)
            disp([' * Preprocessing of ' all_folders(ii).folder filesep all_folders(ii).name]);
            preprocessSession('basepath',pwd,'analogChannelsList',analogChannelsList,'spikeSort',spikeSort,'getPos',getPos, 'cleanArtifacts',cleanArtifacts,...
                'medianSubstr',medianSubstr,'tracking_pixel_cm',tracking_pixel_cm,'sessionSummary',sessionSummary,'digitalChannelsList',digitalChannelsList);
        else 
            disp(['Spikiping ' all_folders(ii).folder filesep all_folders(ii).name]);
        end
        
    end
end


end