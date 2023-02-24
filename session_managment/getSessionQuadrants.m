function [quadrants] = getSessionQuadrants(varargin)
% Computes the quadrants that are on and the HM response (good or bad
% trial)
% INPUTS
%
%
%
%
%
%
%
%
% OUTPUT
%   'quadrants'
%
% Pablo Abad 2022.
%
% ==== POSITIONS ====
% 1: Up left quadrant
% 2: Bottom right quadrant
% 3: Bottom left quadrant
% 4: Up right quadrant
%
% ==== ANSWERS =====
% 7: Up left quadrant
% 3: Bottom right quadrant
% 1: Bottom left quadrant 
% 9: Up right quadrant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse options
p = inputParser;

addParameter(p,'basepath',pwd,@istruct);
addParameter(p,'sr',30000,@isnumeric);
addParameter(p,'force',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
sr = p.Results.sr;
force = p.Results.force;

cd(basepath)
basename = basenameFromBasepath(basepath);

if ~isempty([basenameFromBasepath(basepath) '.Quadrants.Behavior.mat']) && ~force
    file = dir([basenameFromBasepath(basepath) '.Quadrants.Behavior.mat']);
    load(file.name);
end

if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
    load(strcat(basename,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:length(MergePoints.foldernames)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*data_quadrants.mat']))
           cd([basepath filesep MergePoints.foldernames{ii}]);
           fprintf('Computing quadrants in %s folder \n',MergePoints.foldernames{ii});
           tempQuadrants{count} = getQuadrants();
           quadrantsFolder(count) = ii;
           count = count + 1;
        end
    end
end

cd(basepath);
%% Concatenate and sync timestamps
quadrants = []; ts = []; subSessions = []; maskSessions = []; answer = []; position = []; choice = [];
if count > 1
    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(quadrantsFolder)
            if strcmpi(MergePoints.foldernames{quadrantsFolder(ii)},tempQuadrants{ii}.folder)
                sumTs = tempQuadrants{ii}.ts + MergePoints.timestamps(quadrantsFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(quadrantsFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
                
                answer = [answer; tempQuadrants{ii}.answer];
                position = [position; tempQuadrants{ii}.position];
                choice = [choice; tempQuadrants{ii}.choice];
                
                
            else
                error('Folders name does not match!!');
            end
        end
    else
        warning('No MergePoints file found. Concatenating timestamps...');
        for ii = 1:length(quadrantsFolder)
            sumTs = max(ts)+ tempQuadrants{ii}.ts;
            subSessions = [subSessions; [sumTs(1) sumTs(end)]];
            ts = [ts; sumTs];
        end
    end
    
    
    % Concatenating quadrants fields ...
    
    for ii = 1:size(tempQuadrants,2)
    	quadrants.folder{ii} = tempQuadrants{ii}.folder;
        quadrants.performance{ii} = tempQuadrants{ii}.performance;
    end
    
end

% SAVE OUTPUT
quadrants.ts = ts';
quadrants.answer = answer;
quadrants.position = position;
quadrants.choice = choice;

save([basepath filesep basenameFromBasepath(basepath) '.Quadrants.Behavior.mat'],'quadrants');

end