function [fiber] = getSessionFiberPhotometry(varargin)
% [fiber] = getSessionFiberPhotometry()
%
% Get data from a fiber phtometry experiments (.doric file)
% 
% INPUTS
%
%   basepath 
%   green: to load gren signal 
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
if ~isempty(dir([basepath filesep '*fiberPhotometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*fiberPhotometry.mat']); 
    load(file.name)
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
            tempFiber{count} = getFiberPhotometry();
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

    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(fiberFolder)
            if strcmpi(MergePoints.foldernames{fiberFolder(ii)},tempFiber{ii}.folder)
                if isfield(tempFiber{ii},'red_1')
    
                    sumTs = tempFiber{ii}.red_1.timestamps + MergePoints.timestamps(fiberFolder(ii),1);
                    subSessions = [subSessions; MergePoints.timestamps(fiberFolder(ii),1:2)];
                    maskSessions = [maskSessions; ones(size(sumTs))*ii];
                    sumOriginal_ts = tempFiber{ii}.red_1.fiber_timestamps + MergePoints.timestamps(fiberFolder(ii),1);
    
                    ts = [ts; sumTs];
                    original_ts = [original_ts; sumOriginal_ts];
                end
            else
                error('Folders name do not match!!');
            end
        end
    else
        error('No MergePoints file found...');
    end

    % Concatenating fiber fields
    
    green_ts = []; red_ts = []; isosbestic_ts = [];
    green_data = []; red_data = []; isosbestic_data = [];
    red_1_data = []; red_2_data = [];
    red_1_AF_F = []; red_2_AF_F = []; green_AF_F = []; 
    greenL_AF_F = []; greenR_AF_F = [];
    green_zscore = [];
    fiber_sr = [];
    folder = [];
    
    % for ii = 1:size(tempFiber,2)
    %     green_data = [green_ts; tempFiber{ii}.green.data];
    %     red_data = [red_data; tempFiber{ii}.red.data];
    %     isosbestic_data = [isosbestic_data; tempFiber{ii}.isosbestic.data];
    % 
    %     green_AF_F = [green_AF_F; tempFiber{ii}.green.AF_F];
    %     green_zscore = [green_zscore; tempFiber{ii}.green.zscore];
    %     fiber_sr = [fiber_sr; tempFiber{ii}.green.sr];
    %     folder{ii} = tempFiber{ii}.folder;
    % end
    % 
    % fiber.green.data = green_data;
    % fiber.green.AF_F = green_AF_F;
    % fiber.green.zscore = green_zscore;
    % 
    % fiber.red.data = red_data;
    % fiber.isosbestic.data = isosbestic_data;
    % 
    % fiber.timestamps = ts;
    % fiber.sr = fiber_sr(1);
    % fiber.folder = folder;

    for ii = 1:size(tempFiber,2)
        if isfield(tempFiber{ii},'red_1')
            red_1_AF_F = [red_1_AF_F; tempFiber{ii}.red_1.AF_F];
            red_2_AF_F = [red_2_AF_F; tempFiber{ii}.red_2.AF_F];
            % green_AF_F = [green_AF_F; tempFiber{ii}.green.AF_F];
            greenL_AF_F = [greenL_AF_F; tempFiber{ii}.greenL.AF_F];
            greenR_AF_F = [greenR_AF_F; tempFiber{ii}.greenR.AF_F];
            fiber_sr(ii) = [tempFiber{ii}.sr];
            folder{ii} = tempFiber{ii}.folder;
        end
  
    end
    
    fiber.red_1.AF_F = red_1_AF_F;
    fiber.red_2.AF_F = red_2_AF_F;    
    fiber.greenL.AF_F = greenL_AF_F;
    fiber.greenR.AF_F = greenR_AF_F;

    fiber.timestamps = ts;
    fiber.fiber_timestamps = original_ts;
    fiber.sr = fiber_sr(1);
    fiber.folder = folder;
    
    if saveMat
        save([basepath filesep basenameFromBasepath(basepath) '.FiberPhotometry.mat'],'fiber');
    end

else
    warning('No fiber photometry data available!');
    fiber = [];
end

% Save as a CSV file

end