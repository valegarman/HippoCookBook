function [stim] = getSessionStimulation(varargin)
%
%
%
%
%
%
%
%
%
%
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*stimulation.events.mat'])) || forceReload
    disp('stimulation already detected! Loading file.');
    file = dir([basepath filesep '*stimulation.events.mat']);
    load(file.name);
    return
end


cd(basepath);
basename = basenameFromBasepath(basepath);

if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
    load(strcat(basename,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.xml*']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            tempStim{count} = getStimulation();
            stimFolder{count} = ii;
            count = count +1;
        end
    end
    cd(basepath);
else
    error('Missing MergePoints, quitting...');
end

%% Concatenate and sync timestamps
if count > 1 % if traking
    ts = []; subSessions = []; maskSessions = [];
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(stimFolder)
            if strcmpi(MergePoints.foldernames{stimFolder{ii}},tempStim{ii}.folder)
                sumTs = tempStim{ii}.ts + MergePoints.timestamps(stimFolder{ii},1);
                subSessions = [subSessions; MergePoints.timestamps(stimFolder{ii},1:2)];
                maskSessions = [maskSessions ones(size(sumTs))*ii];
                ts = [ts sumTs];
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
    values = [];
    for ii = 1:size(tempStim,2) 
        values = [values tempStim{ii}.values];
        folder{ii} = tempStim{ii}.folder; 
    end

    stim.ts = ts;
    stim.subSessions = subSessions;
    stim.maskSessions = maskSessions;
    stim.values = values;
    stim.folders = folder;
    
    %% save 
    if saveMat
        save([basepath filesep basenameFromBasepath(basepath) '.stimulation.events.mat'],'stim');
    end
else
    warning('No tracking available!');
    tracking = [];
end

end
