
function [tracking] = getSessionTracking(varargin)
%
% Gets position trackign for each sub-session and concatenate all of them so they are 
% aligned with LFP and spikes. Default is recording with Basler, and requiere avi videos 
% and at least one tracking LED. There is an alternative in case OptiTrack was used. 
% Needs to be run in main session folder. 
%
% USAGE
%
%   [tracking] = getSessionTracking(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   roiLED         - 2 x R. 'manual' for drawing the ROI.
%   roisPath       - provide a path with ROI mat files ('roiTRacking.mat'
%                   and 'roiLED.mat'). By default try to find it in
%                   basePath or upper folder.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveMat        - default true
%   forceReload    - default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   position.y               - y position in cm/ normalize
%   timestamps      - in seconds, if Basler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%       only for OptiTrack
%   position.z 
%   orientation.x
%   orientation.y
%   orientation.z

%   HISTORY:
%     - Manuel Valero 2019
%     - Added OptiTrack support: 5/20, AntonioFR (STILL NEEDS TESTING)
%     - Added AnyMaze support: 3/22 Pablo Abad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'convFact',0.1149,@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'roisPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',true,@islogical);
addParameter(p,'LED_threshold',0.98,@isnumeric);
addParameter(p,'tracking_ttl_channel',[],@isnumeric);
addParameter(p,'leftTTL_reward',[],@isnumeric);
addParameter(p,'rightTTL_reward',[],@isnumeric);
addParameter(p,'homeTtl',[],@isnumeric);
addParameter(p,'tracking_software','dlc'); % Options are: 'basler', 'anymaze', 'dlc' (default).
addParameter(p,'interpolate_misstrackings',true);

parse(p,varargin{:});
basepath = p.Results.basepath;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
roisPath = p.Results.roisPath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
LED_threshold = p.Results.LED_threshold;
tracking_ttl_channel = p.Results.tracking_ttl_channel;
leftTTL_reward = p.Results.leftTTL_reward;
rightTTL_reward = p.Results.rightTTL_reward;
homeTtl = p.Results.homeTtl;
tracking_software = p.Results.tracking_software;
interpolate_misstrackings = p.Results.interpolate_misstrackings;



%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

if strcmpi(tracking_software,'basler')
    %% Basler tracking
    cd(basepath); cd ..; upBasepath = pwd; cd(basepath);
    if isempty(roisPath)
        if exist([basepath filesep 'roiTracking.mat'],'file') || ...
            exist([basepath filesep 'roiLED.mat'],'file')
                roisPath = basepath;
                try load([roisPath filesep 'roiLED.mat'],'roiLED'); end
                load([roisPath filesep 'roiTracking.mat'],'roiTracking');
        elseif exist([upBasepath filesep 'roiTracking.mat'],'file') || ...
            exist([upBasepath filesep 'roiLED.mat'],'file')
                roisPath = upBasepath;
                try load([roisPath filesep 'roiLED.mat'],'roiLED'); end
                load([roisPath filesep 'roiTracking.mat'],'roiTracking');
        end   
    end

    %% Find subfolder recordings
    cd(basepath);
    basename = basenameFromBasepath(basepath);
    %C = strsplit(sessionInfo.session.name,'_');
    %sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
             if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*Basler*avi']))   
                cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
                 tempTracking{count}= LED2Tracking([],'convFact',convFact,'roiTracking',...
                     roiTracking,'roiLED',roiLED,'forceReload',forceReload,'thresh',LED_threshold,'basler_ttl_channel',tracking_ttl_channel,'leftTTL_reward',leftTTL_reward,'rightTTL_reward',rightTTL_reward); % computing trajectory
                 trackFolder(count) = ii; 
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('missing MergePoints, quiting...');
    end
elseif strcmpi(tracking_software,'anymaze')
    % Find subfolder recordings
    cd(basepath);
    basename = basenameFromBasepath(basepath);
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*_filtered.csv*']))
                cd([basepath filesep MergePoints.foldernames{ii}]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
                tempTracking{count} = anyMazeTracking([],[]);
                trackFolder(count) = ii;
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('Missing MergePoints, quitting...');
    end

elseif strcmpi(tracking_software,'dlc')
    % Deep lab cut tracking
    cd(basepath);
    basename = basenameFromBasepath(basepath);
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*tracking*_crop.avi']))
                cd([basepath filesep MergePoints.foldernames{ii}]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
                 tempTracking{count} = dlc_tracking('leftTtl_reward',leftTTL_reward,'rightTtl_reward',rightTTL_reward,'homeTtl',homeTtl,'forceReload',forceReload,'interpolate_misstrackings',interpolate_misstrackings);
                trackFolder(count) = ii;
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('Missing MergePoints, quitting...');
    end

end
%% Concatenate and sync timestamps
if count > 1 % if traking
    ts = []; subSessions = []; maskSessions = []; originalts = [];
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(trackFolder)
            if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
                sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
                if strcmpi(tracking_software,'anymaze')
                    sumOriginalTs = tempTracking{ii}.originalTimestamps + MergePoints.timestamps(trackFolder(ii),1);
                    originalts = [originalts; sumOriginalTs];
                end
            else
                error('Folders name does not match!!');
            end
        end
    else
        warning('No MergePoints file found. Concatenating timestamps...');
        for ii = 1:length(trackFolder)
            sumTs = max(ts)+ tempTracking{ii}.timestamps;
            subSessions = [subSessions; [sumTs(1) sumTs(end)]];
            ts = [ts; sumTs];
        end
    end

    % Concatenating tracking fields...
    x = []; y = []; likelihood = [];
    x_back = []; y_back = []; likelihood_back = [];
    x_tail1 = []; y_tail1 = []; likelihood_tail1 = [];
    x_tail2 = []; y_tail2 = []; likelihood_tail2 = [];
    
    
    folder = []; samplingRate = []; description = [];
    velocity = []; acceleration = []; 
    avFrame = []; apparatus = []; zone = []; pixelsmetre = [];
    for ii = 1:size(tempTracking,2) 
        x = [x; tempTracking{ii}.position.x]; 
        y = [y; tempTracking{ii}.position.y];
        likelihood = [likelihood; tempTracking{ii}.position.likelihood];
        x_back = [x_back; tempTracking{ii}.position_back.x];
        y_back = [y_back; tempTracking{ii}.position_back.y];
        likelihood_back = [likelihood_back; tempTracking{ii}.position_back.likelihood];
        x_tail1 = [x_tail1; tempTracking{ii}.position_tail1.x];
        y_tail1 = [y_tail1; tempTracking{ii}.position_tail1.y];
        likelihood_tail1 = [likelihood_tail1; tempTracking{ii}.position_tail1.likelihood];
        x_tail2 = [x_tail2; tempTracking{ii}.position_tail2.x];
        y_tail2 = [y_tail2; tempTracking{ii}.position_tail2.y];
        likelihood_tail2 = [likelihood_tail2; tempTracking{ii}.position_tail2.likelihood];


            
        velocity = [velocity; tempTracking{ii}.velocity];
        acceleration = [acceleration; tempTracking{ii}.acceleration];
        folder{ii} = tempTracking{ii}.folder; 
        samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
        description{ii} = tempTracking{ii}.description;
        if strcmpi(tracking_software,'anymaze')
            avFrame{ii} = tempTracking{ii}.avFrame;
            apparatus{ii} = tempTracking{ii}.apparatus;
            if isfield(tempTracking{ii},zone)
                zone{ii} = tempTracking{ii}.zone;
            else
                zone{ii} = [];
            end
            pixelsmetre{ii} = tempTracking{ii}.pixelsmetre;
        end
    end

    tracking.position.x = x;
    tracking.position.y = y;
    tracking.folders = folder;
    tracking.samplingRate = samplingRate;
    tracking.timestamps = ts;
    tracking.events.subSessions =  subSessions;
    tracking.events.subSessionsMask = maskSessions;
    % Create speed structure .mat
    
    speed.velocity = velocity;
    speed.acceleration = acceleration;
    speed.timestamps = ts;
    speed.folders = folder;
    
    if strcmpi(tracking_software,'anymaze')
        tracking.avFrame = avFrame;
        tracking.apparatus = apparatus;
        tracking.zone = zone;
        tracking.pixelsmetre = pixelsmetre;
    end


    %% save tracking 
    if saveMat
        save([basepath filesep basenameFromBasepath(basepath) '.Tracking.Behavior.mat'],'tracking');
        save([basepath filesep basenameFromBasepath(basepath), '.Speed.events.mat'],'speed');
    end
else
    warning('No tracking available!');
    tracking = [];
end

end

