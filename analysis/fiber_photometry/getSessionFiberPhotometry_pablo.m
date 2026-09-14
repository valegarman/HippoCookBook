function [fiber] = getSessionFiberPhotometry_pablo(varargin)
% [fiber] = getSessionFiberPhotometry()
%
% Get data from a fiber phtometry experiments (.doric file)
% 
% INPUTS
%
%   basepath 
%   green: to load green signal 
%   red: to load red signal
%   isobestic: to load isobestic signal
%   forceReload
%   saveMat     - default, true

% Inputs
p = inputParser();

addParameter(p,'basepath',pwd);
addParameter(p,'saveMat',true);
addParameter(p,'force',false);

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
force = p.Results.force;

%% In case already exists
if ~isempty(dir([basepath filesep '*FiberPhotometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*FiberPhotometry.mat']); 
    load(file.name)
    return;
end


% Load session
session = loadSession();

% Find subfolders recordings
cd(basepath);
basename = basenameFromBasepath(basepath);
if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
    load(strcat(basename,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.doric']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
            fprintf('Computing fiber photometry in %s folder \n',MergePoints.foldernames{ii});
            % tempFiber{count} = getFiberPhotometry();
            % tempFiber{count} = getFiberPhotometry_temp_Nacho('force',true);
            % tempFiber{count} = getFiberPhotometry_temp_Nacho_v2('force',true);
            tempFiber{count} = getFiberPhotometry_pablo('force',true);

            fiberFolder(count) = ii;
            count = count +1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end

%% Concatenate and sync timestamps

if count > 1 % if fiber recording

    ts = []; subSessions = []; maskSessions = []; original_ts = [];
    % Event Features
    startpoint_timestamp = []; endpoint_timestamp = []; peak_timestamp = [];
    timestamps = []; aligned_timestamps = [];

    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(fiberFolder)
            if strcmpi(MergePoints.foldernames{fiberFolder(ii)},tempFiber{ii}.folder)
                if isfield(tempFiber{ii},'timestamps')
    
                    sumTs = tempFiber{ii}.timestamps + MergePoints.timestamps(fiberFolder(ii),1);
                    subSessions = [subSessions; MergePoints.timestamps(fiberFolder(ii),1:2)];
                    maskSessions = [maskSessions; ones(size(sumTs))*ii];
                    sumOriginal_ts = tempFiber{ii}.original_timestamps + MergePoints.timestamps(fiberFolder(ii),1);
    
                    ts = [ts; sumTs];
                    original_ts = [original_ts; sumOriginal_ts];

                    % Event Features
                    sum_startpoint_ts = tempFiber{ii}.eventFeatures.Startpoint_Timestamp + MergePoints.timestamps(fiberFolder(ii),1);
                    startpoint_timestamp = [startpoint_timestamp; sum_startpoint_ts];

                    sum_endpoint_ts = tempFiber{ii}.eventFeatures.Endpoint_Timestamp + MergePoints.timestamps(fiberFolder(ii),1);
                    endpoint_timestamp = [endpoint_timestamp; sum_endpoint_ts];

                    sum_peak_ts = tempFiber{ii}.eventFeatures.PeakTimestamps + MergePoints.timestamps(fiberFolder(ii),1);
                    peak_timestamp = [peak_timestamp; sum_peak_ts];

                end
            else
                error('Folders name do not match!!');
            end
        end
    else
        error('No MergePoints file found...');
    end

    % Concatenating fiber fields

    timestamps = []; original_timestamps = []; isosbestic = []; green = []; red = []; 
    isosbestic_original = []; green_original = []; red_original = [];
    % Event features
    event_id = []; event_frequency = [];  dff_smoothed = []; aligned_dff_smoothed = [];  duration = []; amplitude = []; fluorescent_change = [];
    time_to_peak = []; auc = [];


    
    for ii = 1:size(tempFiber,2)

        sr{ii} = tempFiber{ii}.sr;
        folder{ii} = tempFiber{ii}.folder;  

        if isfield(tempFiber{ii},'red')          
            red = [red; tempFiber{ii}.red];
            red_original = [red_original; tempFiber{ii}.red_original.data];

        end

        if isfield(tempFiber{ii},'green')
   
            green = [green; tempFiber{ii}.green];
            green_original = [green_original; tempFiber{ii}.green_original.data];

        end

        if isfield(tempFiber{ii},'iso')

            isosbestic = [isosbestic; tempFiber{ii}.iso];
            isosbestic_original = [isosbestic_original; tempFiber{ii}.iso_original.data];
        end

        % Preprocessing
        preprocessing{ii} = tempFiber{ii}.preprocessing;

        % eventFeatures
        event_id = [event_id; tempFiber{ii}.eventFeatures.Event_ID];
        event_frequency = [event_frequency; tempFiber{ii}.eventFeatures.Event_Frequency]; 
        
        dff_smoothed = [dff_smoothed; tempFiber{ii}.eventFeatures.DFF_Smoothed]; 
        aligned_dff_smoothed = [aligned_dff_smoothed; tempFiber{ii}.eventFeatures.Aligned_DFF_Smoothed]; 

        duration = [duration; tempFiber{ii}.eventFeatures.Duration]; 
        amplitude = [amplitude; tempFiber{ii}.eventFeatures.Amplitude]; 
        fluorescent_change = [fluorescent_change; tempFiber{ii}.eventFeatures.FluorescenceChange];
        time_to_peak = [time_to_peak; tempFiber{ii}.eventFeatures.TimeToPeak]; 
        auc = [auc; tempFiber{ii}.eventFeatures.AUC];

    end

    fiber = [];

    fiber.timestamps = ts;
    fiber.original_timestamps = original_ts;
    fiber.red = red;
    fiber.red_original = red_original;
    fiber.green = green;
    fiber.green_original = green_original;
    fiber.iso = isosbestic;
    fiber.iso_original = isosbestic_original;

    fiber.preprocessing = preprocessing;
    % Event Features
    fiber.eventFeatures.startpoint_timestamp = startpoint_timestamp;
    fiber.eventFeatures.endpoint_timestamp = endpoint_timestamp;
    fiber.eventFeatures.peak_timestamp = peak_timestamp; 
    fiber.eventFeatures.event_id = event_id; 
    fiber.eventFeatures.event_frequency = event_frequency;
    fiber.eventFeatures.dff_smoothed = dff_smoothed;
    fiber.eventFeatures.aligned_dff_smoothed = aligned_dff_smoothed;
    fiber.eventFeatures.duration = duration;
    fiber.eventFeatures.amplitude = amplitude;
    fiber.eventFeatures.fluorescent_change = fluorescent_change;
    fiber.eventFeatures.time_to_peak = time_to_peak;
    fiber.eventFeatures.auc = auc;

    fiber.sr = sr{1};
    fiber.folder = folder;

    fiber.events.subSessions = subSessions;
    fiber.events.maskSessions = maskSessions;
    
    if saveMat
        save([basepath filesep basenameFromBasepath(basepath) '.FiberPhotometry.mat'],'fiber');
    end

else
    warning('No fiber photometry data available!');
    fiber = [];
end

close all;

% Save as a CSV file

end