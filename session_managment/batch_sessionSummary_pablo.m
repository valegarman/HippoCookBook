function  batch_sessionSummary_pablo(varargin)

%   Runs a list of preliminary descriptive analysis throughout all sessions 
%   without a SummaryFigures folder on them.
% Runs a list of preliminary descriptive analysis 
%
% USAGE 
%   computeSessionSummary(varargin)
%
% INPUT (optional)
%   basepath         By default pwd
%   listOfAnalysis   Summary analysis that will be computed (as a cell with
%                        strings). Default 'all'. Possible values: 'spikes',
%                        'analogPulses', 'digitalPulses',
%                        'downStates','ripples', 'tMazeBehaviour',
%                        'linearMazeBehaviour', 'thetaModulation'
%   exclude         Cell with strings with the list of analysis to exclude
%  
%   (specific analysis optionss)
%   excludeShanks           Default []
%   analogChannelsList      Array of channel to perform 'analogPulses' psth and csd. Default 'all'
%   digitalChannelsList     Array of channel to perform 'digitalPulses' psth and csd. Default 'all'
%
% Manu Valero-BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'analysisPath',[],@isdir);
addParameter(p,'listOfAnalysis','all',@iscellstr);
addParameter(p,'exclude',[],@iscellstr);
addParameter(p,'excludeShanks',[],@isnumeric);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'force',false,@islogical);
addParameter(p,'anyMaze',true,@islogical);
addParameter(p,'pixelsPerCm',2.5,@isnumeric);


parse(p,varargin{:});
basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;
tracking_pixel_cm = p.Results.tracking_pixel_cm;
force = p.Results.force;
anyMaze = p.Results.anyMaze;
pixelsPerCm = p.Results.pixelsPerCm;

% batch processing...
all_folders = dir(basepath);
for ii = 1:size(all_folders,1)
    if all_folders(ii).isdir && ~strcmpi(all_folders(ii).name,'.') && ~strcmpi(all_folders(ii).name,'..')
        cd([all_folders(ii).folder filesep all_folders(ii).name]);
        summaryFolder = dir('*SummaryFigures*');
        if isempty(summaryFolder) || force
            disp(['Computing sessionSummary of ' all_folders(ii).folder filesep all_folders(ii).name]);
            computeSessionSummary_pablo('basepath',pwd,'listOfAnalysis',listOfAnalysis,'exclude',exclude,'excludeShanks',excludeShanks, 'analogChannelsList',analogChannelsList,...
                'digitalChannelsList',digitalChannelsList,'tracking_pixel_cm',tracking_pixel_cm,'anyMaze',anyMaze,'pixelsPerCm',pixelsPerCm);
        else 
            disp(['Spikiping ' all_folders(ii).folder filesep all_folders(ii).name]);
        end
        
    end
end

end