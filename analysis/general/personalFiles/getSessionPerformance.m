function [performance] = getSessionPerformance(varargin)
%
% Gets session performance and variables for different behavioral paradigms
%
% USAGE
%
%   [performance] = getSessionPerformance(varargin)
%
% INPUTS
%
% OUTPUT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'includeIntervals',[],@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
includeIntervals = p.Results.includeIntervals;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Performance.Behavior.mat'])) || forceReload
    disp('Performance already detected! Loading file.');
    file = dir([basepath filesep '*Performance.Behavior.mat']);
    load(file.name);
    return
end

% Load tracking
tracking = getSessionTracking();


% Load behavior
behavior = getSessionBehavior();

% Load speed
if ~isempty(dir('*speed.events.mat'))
    file = dir('*speed.events.mat');
    load(file.name);
end

[status] = InIntervals(tracking.timestamps,includeIntervals);
subsessions = unique(tracking.events.subSessionsMask(status));

for ii = 1:length(subsessions)
    
    description = behavior.description{subsessions(ii)};
    
    xPos = tracking.position.x(tracking.events.subSessionsMask == subsessions(ii));
    yPos = tracking.position.y(tracking.events.subSessionsMask == subsessions(ii));

    for jj = 2:length(xPos)-1
        dist(jj) = sqrt((xPos(jj)-xPos(jj-1))^2 + (yPos(jj)-yPos(jj-1))^2);
    end
    
    distance = max(cumsum(dist))/100;
    ts = tracking.timestamps(status);
    start = ts(1);
    stop = ts(end);
    distance_perc = distance / ((stop-start)/60);
    meanSpeed = mean(speed.velocity(tracking.events.subSessionsMask == subsessions(ii)));
    
    if strcmpi(behavior.description{subsessions(ii)},'Open Field')
        performance.OpenField.meanSpeed = meanSpeed;
        performance.OpenField.distance = distance;
        performance.OpenField.distance_perc = distance_perc;
    elseif strcmpi(behavior.description{subsessions(ii)},'YMaze Apparatus')
        performance.YMaze.meanSpeed = meanSpeed;
        performance.YMaze.distance = distance;
        performance.YMaze.distance_perc = distance_perc;
    end
    figure;
    plot(xPos,yPos);
end
