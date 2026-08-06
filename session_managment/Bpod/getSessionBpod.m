function [tracking] = getSessionBpod(varargin)
%
% Gets Bpod (Sanworks) information (events and states) for each sub-session and concatenates.  
%
% USAGE
%
%   [Bpod] = getSessionBpod(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file
%
% OUTPUT
%   Bpod.behaviour output structure, with fields:
%       - performance.
%
%   HISTORY:
%     - Pablo Abad 2026 (NeuralComputationLab).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Bpod.Behavior.mat'])) && forceReload
    disp('Bpod already detected! Loading file.');
    file = dir([basepath filesep '*Bpod.Behavior.mat']);
    load(file.name);
    return
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
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*CentralPokeReward*mat']))   
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
    files = dir('*_2*');
    count = 1;
    for ii = 1:size(files,1)
        if ~isempty(dir([basepath filesep files(ii).name filesep '*CentralPokeReward*mat']))
            cd([basepath filesep files(ii).name]);
            fprintf('Computing Bpod in %s folder \n',files(ii).name);
            tempBpod{count}= getBpod(); 
            BpodFolder(count) = ii; 
            count = count + 1;
        end
    end
    warning('missing MergePoints, quiting...');
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
                maskSessions = [maskSessions; ones(size(sumTs))'*ii];
                ts = [ts; sumTs'];
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