function [tracking,behaviour] = checkTracking(varargin)
%
% Checks if there are conditions where the duration of the behaviour is
% much bigger than other ones.
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
%   folder          
% Pablo Abad, 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'minTime',5,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
minTime = p.Results.minTime;

%% In case tracking already exists 
tracking = getSessionTracking();
behaviour = getSessionBehaviour();

for i=1:length(tracking.apparatus)
    apparatus_name{i} = tracking.apparatus{i}.name;
end
unique_apparatus_name = unique(apparatus_name);

for ii = 1:length(unique_apparatus_name)
    index = find(ismember(apparatus_name,unique_apparatus_name{ii}));
    
    for jj = 1:length(index)
        timestamps(jj) = length(find(tracking.events.subSessionsMask == index(jj)));
    end
    difTime = abs(diff(timestamps))/30/60;
    if difTime > minTime
        % Correcting timestamps...
        [mTime,posmTime] = min(timestamps);
        [maxTime,posMaxTime] = max(timestamps);
        fprintf('Correcting timestamps of %s folder \n', tracking.folders{index(posMaxTime)});
        indexToExclude = zeros(length(tracking.events.subSessionsMask),1);
        indexToExclude(mTime:find(tracking.events.subSessionsMask == index(posMaxTime),1,'last'),1) = 1;
        
        tracking.position.x(find(indexToExclude)) = [];
        tracking.position.y(find(indexToExclude)) = [];
        tracking.timestamps(find(indexToExclude)) = [];
        tracking.events.subSessionsMask(find(indexToExclude)) = [];
        
        behaviour.timestamps(find(indexToExclude)) = [];
        behaviour.position.x(find(indexToExclude)) = [];
        behaviour.position.y(find(indexToExclude)) = [];
        behaviour.events{index(posMaxTime)}.subSessionMask = tracking.events.subSessionsMask(tracking.events.subSessionsMask == index(posMaxTime));
        behaviour.maps{index(posMaxTime)} = behaviour.maps{index(posMaxTime)}(1:mTime,:);
    end
end




%% save tracking 
if saveMat
    save([basepath filesep basenameFromBasepath(basepath) '.Tracking.Behavior.mat'],'tracking');
    save([basepath filesep basenameFromBasepath(basepath) '.Behavior.mat'],'behaviour');
end

end
