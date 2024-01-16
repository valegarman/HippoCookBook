function [tracking] = getSessionTracking_SR(varargin)
%
% Gets position trackign for simultaneous recordings and save it in respective folder 
%
% USAGE
%
%   [tracking] = getSessionTracking_SR(varargin)
%
% INPUTS
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
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'anyMaze',false,@islogical);
addParameter(p,'LED_threshold',0.98,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
roisPath = p.Results.roisPath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
anyMaze = p.Results.anyMaze;
LED_threshold = p.Results.LED_threshold;


% Find subfolder recordings
cd(basepath);
basename = basenameFromBasepath(basepath);
if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
    load(strcat(basename,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.csv*']))
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



end


