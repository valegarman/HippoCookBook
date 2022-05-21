

function findOverlapStimInterval(varargin)

p = inputParser();
addParameter(p,'basepath',pwd);
addParameter(p,'optoStim_names',pwd);

parse(p,varargin{:});
basepath = p.Results.basepath;
optoStim_names = p.Results.optoStim_names;

% Winnie Yang, 2022
% NOTICE, only work if there are two optical fibers for now
%%
basename = basenameFromBasepath(basepath);

%optoStim_names = {'digitalIn1', 'digitalIn2'};
%%



optoStim_events = {};
for nn = 1:length(optoStim_names)
    load([basename,'.',optoStim_names{nn}, '.events.mat']);
    optoStim_events{nn} = events;
end



%%

OVERLAP_1 = [];
OVERLAP_2 = [];

% START TIMESTAMP
[overlapEvent_start1,overlapEvent_start2,~] = InIntervals(optoStim_events{1}.timestamps(:,1), optoStim_events{2}.timestamps);
overlapEvent_start2 = overlapEvent_start2(overlapEvent_start2>0);
overlapEvent_start1 = find(overlapEvent_start1>0);

%STOP TIMESTAMP
[overlapEvent_stop1,overlapEvent_stop2,~] = InIntervals(optoStim_events{1}.timestamps(:,2), optoStim_events{2}.timestamps);
overlapEvent_stop2 = overlapEvent_stop2(overlapEvent_stop2>0);
overlapEvent_stop1 = find(overlapEvent_stop1>0);

overlapEvents_1 = union(overlapEvent_start1,overlapEvent_stop1);
overlapEvents_2 = union(overlapEvent_start2,overlapEvent_stop2);

OVERLAP_1 = [OVERLAP_1,overlapEvents_1];
OVERLAP_2 = [OVERLAP_2,overlapEvents_2];


%% REPEAT, OTHER DIRECTION 

%START TIMESTAMP
[overlapEvent_start2,overlapEvent_start1,~] = InIntervals(optoStim_events{2}.timestamps(:,1), optoStim_events{1}.timestamps);
overlapEvent_start1 = overlapEvent_start1(overlapEvent_start1>0);
overlapEvent_start2 = find(overlapEvent_start2>0);

%STOP TIMESTAMP
[overlapEvent_stop2,overlapEvent_stop1,~] = InIntervals(optoStim_events{2}.timestamps(:,2), optoStim_events{1}.timestamps);
overlapEvent_stop1 = overlapEvent_stop1(overlapEvent_stop1>0);
overlapEvent_stop2 = find(overlapEvent_stop2>0);

overlapEvents_1 = union(overlapEvent_start1,overlapEvent_stop1);
overlapEvents_2 = union(overlapEvent_start2,overlapEvent_stop2);

OVERLAP_1 = [OVERLAP_1,overlapEvents_1];
OVERLAP_2 = [OVERLAP_2,overlapEvents_2];


events_1 = 1:length(optoStim_events{1}.timestamps);
events_2 = 1:length(optoStim_events{2}.timestamps);


events_1(OVERLAP_1) = 0;
noOverlap_1 = find(events_1>0);
overlap_1 = find(events_1==0);

events_2(OVERLAP_2) = 0;
noOverlap_2 = find(events_2>0);
overlap_2 = find(events_2==0);




events_1 ={};
events_2 ={};
events_1.noOverlap = noOverlap_1';
events_1.overlap = overlap_1';

events_2.noOverlap = noOverlap_2';
events_2.overlap = overlap_2';
% events_1.timestamps = optoStim_events{1}.timestamps(noOverlap_1,:);
% events_2.timestamps = optoStim_events{2}.timestamps(noOverlap_2,:);
% 
% 
% events_1.duration = optoStim_events{1}.duration(noOverlap_1)';
% events_2.duration = optoStim_events{2}.duration(noOverlap_2)';

save([basename,'.optogenetic_overlapEvents_1.mat'],'events_1');
save([basename,'.optogenetic_overlapEvents_2.mat'],'events_2');









