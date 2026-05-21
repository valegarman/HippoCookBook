function [optoPulses] = getFiberOptogeneticPulses(varargin)

%       [optoPulses] = getFiberOptogeneticPulses(varargin)
%
%   Gets optogenetic pulses and detects blocks of stimulation (to only keep
%   the first TTL).

% INPUTS 
%
%   <OPTIONALS>
%   digitalChannelList          List of digital channels with light pulses. By default, []%
%
%
% OUTPUT
%   optoPulses
%
%
%   Develop by Pablo Abad. Neural Computational Lab 2025


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'force',true);
addParameter(p,'plt',true);
addParameter(p,'savePlot',true);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'digitalChannelsList',NaN);
addParameter(p,'minNumberOfPulses',99,@isnumeric);
addParameter(p,'duration_round_decimal',3,@isscalar);
addParameter(p,'minDuration',0.004,@isnumeric); % 4 ms
addParameter(p,'maxNumberOfPulses',5000,@isnumeric);
addParameter(p,'min_inter_block',30); % 20 seconds between blocks
addParameter(p,'min_inter_repetitions',20);


parse(p,varargin{:})

basepath = p.Results.basepath;
force = p.Results.force;
plt = p.Results.plt;
savePlot = p.Results.savePlot;
restrict_to = p.Results.restrict_to;
digitalChannelsList = p.Results.digitalChannelsList;
minNumberOfPulses = p.Results.minNumberOfPulses;
duration_round_decimal = p.Results.duration_round_decimal;
minDuration = p.Results.minDuration;
maxNumberOfPulses = p.Results.maxNumberOfPulses;
min_inter_block = p.Results.min_inter_block;
min_inter_repetitions = p.Results.min_inter_repetitions;

% Deal with inputs
prevPath = pwd;
cd(basepath);

% Load session
session = loadSession();

if exist([session.general.name '.optogeneticFiberPulses.events.mat']) && ~force
    disp(['optoPulses already computed for', session.general.name,'.Loading file.']);
    load([session.general.name '.optogeneticFiberPulses.events.mat']);
end

%% Digital pulses
   
if isnan(digitalChannelsList)
    try
        session = loadSession(basepath);
        digitalChannelsList = session.analysisTags.digital_optogenetic_channels;
    catch
        warning('There is a problem with the digital channels...');
    end
end

pulsesDigital.timestamps = []; pulsesDigital.digitalChannelsList = [];
if ~isempty(digitalChannelsList)
    digitalIn = getDigitalIn;
    if ~isempty(digitalIn)
        for ii = 1:length(digitalChannelsList)
            pulsesDigital.timestamps = [pulsesDigital.timestamps; digitalIn.ints{digitalChannelsList(ii)}];
            pulsesDigital.digitalChannelsList = [pulsesDigital.digitalChannelsList; [ones(size(digitalIn.ints{digitalChannelsList(ii)},1),1) * digitalChannelsList(ii)]];
        end
    end
end

pulses.timestamps = [pulsesDigital.timestamps];  % combine pulses
if isempty(pulses.timestamps)
    optogeneticResponses = [];
    return 
end

pulses.channel = [pulsesDigital.digitalChannelsList];  % combine pulses
pulses.digitalChannelsList = [pulsesDigital.digitalChannelsList];  % 
pulses.duration = round(pulses.timestamps(:,2) - pulses.timestamps(:,1),3);  % 
pulses.isDigital = [ones(size(pulsesDigital.digitalChannelsList))];

% restrict_pulses
% status = InIntervals(pulses.timestamps(:,1),restrict_ints);
% pulses.timestamps = pulses.timestamps(status,:);
% pulses.channel = pulses.channel(status,:);
% pulses.digitalChannelsList = pulses.digitalChannelsList(status,:);
% pulses.duration = pulses.duration(status,:);
% pulses.isDigital = pulses.isDigital(status,:);
% pulses.restricted_intervals = restrict_ints;

pulseDuration = unique(round(pulses.duration,duration_round_decimal)); % because code only codes for channel, we take minimum duration channel for responses
channels = unique(pulses.channel); % code per channel, channel x duration should be implemented... 
timestamps_recording = min(pulses.timestamps(:,2)):1/1250:max(pulses.timestamps(:,2));
% pulses condition channels x durations
[m,n] = ndgrid(pulseDuration,channels);
conditions = [m(:),n(:)];

for ii = 1:size(conditions,1)
    conditions(ii,3) = length(find(pulses.duration==conditions(ii,1) & pulses.channel == conditions(ii,2)));
end
notEnoughtPulses = conditions(:,3)<minNumberOfPulses;
conditions(notEnoughtPulses,:) = []; % removing groups of pulses with less number of pulses than defined in 'notEnoughtPulses'
conditions(conditions(:,1)==0,:) = []; % removing pulses with duration shorter than decimal round
nConditions = size(conditions,1);

%% Merge durations that are close enough to be the same pulse duration
for i = 1:length(channels)
    indexes{i} = find(conditions(:,2) == channels(i));
    durations = conditions(find(conditions(:,2) == channels(i)),1);
    index{i} = find(diff(abs(durations)) < minDuration);
end

for i = 1:length(index)
    if ~isempty(index{i})
        for j = 1:length(index{i})
             indexToMerge{i}(j) = find(conditions(:,3) == max(conditions(index{i}(j),3), conditions(index{i}(j)+1,3)),1);
        end
    else
        indexToMerge{i} = [];
    end
end

% Modify pulses.durations values
for i = 1:length(channels)
    if ~isempty(index{i})
        ind = indexes{i}(index{i});
        indToMerge = indexes{i}(indexToMerge{i});
        for j = 1:length(ind)
            dur = find(pulses.duration == conditions(ind(j),1) & pulses.channel == conditions(ind(j),2));
            pulses.duration(dur) = conditions(indToMerge(j),1);
        end
    end
end

% Modify conditions
for i = 1:length(channels)
    if ~isempty(index{i})
        ind = indexes{i}(index{i});
        indToMerge = indexes{i}(indexToMerge{i});
        for j = 1:length(ind)
            conditions(indToMerge(j),3) = conditions(indToMerge(j),3) + conditions(ind(j),3);
        end
    end
end

% Delete the indices
indxs = [];
for i = 1:length(channels)
    indxs = [indxs ;indexes{i}(index{i})];
end
conditions(indxs,:) = [];
nConditions = size(conditions,1);

if nConditions == 2
    if abs(conditions(1,1) - conditions(2,1)) < minDuration
        conditions(2,3) = conditions(2,3) + conditions(1,3);
        minPulse = min(conditions(:,1));
        maxPulse = max(conditions(:,1));
        
        pulses.duration(pulses.duration == minPulse) = maxPulse;
        
        pulseDuration = unique(round(pulses.duration,3)); % because code only codes for channel, we take minimum duration channel for responses
        [m,n] = ndgrid(pulseDuration,channels);
        conditions = [m(:),n(:)];
        for ii = 1:size(conditions,1)
            conditions(ii,3) = length(find(pulses.duration==conditions(ii,1)));
        end
        notEnoughtPulses = conditions(:,3)<minNumberOfPulses;
        conditions(notEnoughtPulses,:) = [];
        nConditions = size(conditions,1);
    end
end

toRemove = find(pulses.duration< minDuration);
pulses.timestamps(toRemove,:) = [];
pulses.channel(toRemove) = [];
pulses.digitalChannelsList(toRemove) = [];
pulses.duration(toRemove) = [];
pulses.isDigital(toRemove) = [];
pulses.nConditions = nConditions;


% Try to detect different blocks of stimulation
d = diff(pulses.timestamps(:,1));
repetition_starts_idx = find(d > min_inter_repetitions) +1; % Find high jumps, so start of a new repetition
repetition_starts_idx = [1; repetition_starts_idx(:)]; % Make sure that we include the first TTL
repetition_starts_times = pulses.timestamps(repetition_starts_idx); % Extract real timestamps

% Let's find how many repetitions per block
block_starts_idx = find(d > min_inter_block) +1; % Find high jumps, so start of a new repetition
block_starts_idx = [1; block_starts_idx(:)]; % Make sure that we include the first TTL
block_starts_times = pulses.timestamps(block_starts_idx); % Extract real timestamps

repetitions_per_block = numel(repetition_starts_idx)/numel(block_starts_idx);
fprintf('Detected %d stimulation blocks and %d repetitions per block\n', numel(block_starts_idx), repetitions_per_block);

% Now we group the pulses per block
blocks = cell(length(block_starts_idx),1);

for ii = 1:length(blocks)

    idx_start = (ii-1)*repetitions_per_block + 1;
    idx_end   = idx_start + repetitions_per_block - 1;
    
    blocks{ii} = repetition_starts_idx(idx_start:idx_end);
    blocks{ii} = repetition_starts_times(idx_start:idx_end);
    
end

% find intervals
lag = 120; % max intrerval between pulses, in seconds
stimulationEpochs(1,1) = pulses.timestamps(1,1);
intPeaks =find(diff(pulses.timestamps(:,1))>lag);
for ii = 1:length(intPeaks)
    stimulationEpochs(ii,2) = pulses.timestamps(intPeaks(ii),2);
    stimulationEpochs(ii+1,1) = pulses.timestamps(intPeaks(ii)+1,1);
end
stimulationEpochs(end,2) = pulses.timestamps(end,2);
pulses.stimulationEpochs = stimulationEpochs;

% Write output

optoPulses = [];
optoPulses = pulses;
if length(blocks) > 1
    optoPulses.trains = [];
    optoPulses.trains.timestamps = cell2mat(cellfun(@(x) x(:).', blocks, 'UniformOutput', false));
    optoPulses.trains.durations = sort(unique(diff(repetition_starts_idx)));
    optoPulses.trains.frequency = sort(unique(diff(repetition_starts_idx)));
end

save([session.general.name,'.optogeneticFiberPulses.events.mat'],'optoPulses');

end