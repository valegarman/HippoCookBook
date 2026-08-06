function [fiberOptogeneticResponse] = getFiberOptogeneticResponse_pablo(varargin)

%       [fiberOptogeneticResponse] = getFiberOptogeneticResponse(timestamps,varargin)
%
%   Computes fiber photometry response and statistics in
%   specific timestamps. 

% INPUTS 
% ts: timestamps for computing fiber photometry responses
%
%   <OPTIONALS>
%   digitalChannelList          List of digital channels with light pulses. By default, []
%   fiber                       Fiber photometry structure. If not provided, tries getSessionFiberPhotometry_temp
%   numRep                      For bootstraping, default 500. If 0, no
%       bootstrapping.
%   winSize                     In seconds, default [-5 5]
%   
%   
%   

%
%   event_ints - interval around events timestamps to compute fiber
%       responses. Default: [-1 5]
%   baseline_ints - interval before even timestamps to compute baseline. Default: [-8 -1]

%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'win'  window over which compute fiber photometry activity
%    =========================================================================
%
% OUTPUT


%   Develop by Pablo Abad. Neural Computational Lab 2024
warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'fiber',[]);
addParameter(p,'winSize',[-3 3],@isnumeric);
addParameter(p,'reload_fiber',false);
addParameter(p,'force',true);
addParameter(p,'event_ints',[0 3]); 
addParameter(p,'baseline_ints',[-3 0]);
addParameter(p,'plt',true);
addParameter(p,'savePlot',true);


addParameter(p,'numRep',500);
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
fiber = p.Results.fiber;
winSize = p.Results.winSize;
reload_fiber = p.Results.reload_fiber;
force = p.Results.force;
event_ints = p.Results.event_ints;
baseline_ints = p.Results.baseline_ints;
plt = p.Results.plt;
savePlot = p.Results.savePlot;

numRep = p.Results.numRep;
restrict_to = p.Results.restrict_to;
digitalChannelsList = p.Results.digitalChannelsList;
minNumberOfPulses = p.Results.minNumberOfPulses;
duration_round_decimal = p.Results.duration_round_decimal;
minDuration = p.Results.minDuration;
maxNumberOfPulses = p.Results.maxNumberOfPulses;
min_inter_block = p.Results.min_inter_block;
min_inter_repetitions = p.Results.min_inter_repetitions;




%% Load fiber

if isempty(fiber)
    fiber = getSessionFiberPhotometry_pablo('force',reload_fiber);
end

% Default parameters
warning('Using default parameters!');
win = winSize;
% win_size = round(fiber.sr * win);
event_ints = [0 3];
baseline_ints = [-3 -3+diff(event_ints)]; 
time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);
c_axis = 3;

t_duringPulse = time_vector > event_ints(1) & time_vector < event_ints(2);
t_beforePulse = time_vector > -3 & time_vector < -1;
t_Z = time_vector >= -3 & time_vector < 0;

eventType = 'opto';

% Deal with inputs
prevPath = pwd;
cd(basepath);

% Load session
session = loadSession();

% if exist([session.general.name '.' eventType '_fiber.mat']) && ~force
%     disp(['Fiber already computed for', session.general.name, ' ', eventType, '.Loading file.']);
%     load([session.general.name '.' eventType '_fiber.mat']);
% end

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

%% Fiber analysis

disp('Computing responses...');

if length(blocks) > 1
    disp('Analyzing blocks of stimulation');
    % Analyze fiber
    count = 0;
    for ii = 1:length(blocks)
    
        fprintf('\n Block %3.i / %3.i \n', ii, length(blocks))
    
        pul = blocks{ii};
        if isempty(pul)
            pul = [0];
        end
        if length(pul)>maxNumberOfPulses
            pul = sort(randsample(pul, maxNumberOfPulses));
        end
    
        if isfield(fiber, 'red') 
            red{ii}.responsecurve = [];
            red_corrected{ii}.responsecurve = [];
            red_smooth{ii}.responsecurve = [];

            if fiber.clean_artifacts
                red_clean{ii}.responsecurve = [];
                red_corrected_clean{ii}.responsecurve = [];
                red_smooth_clean{ii}.responsecurve = [];
            end
            disp('Analyzing red');

        end
        if isfield(fiber, 'green') 
            green{ii}.responsecurve = [];
            green_corrected{ii}.responsecurve = [];
            green_smooth{ii}.responsecurve = [];
            if fiber.clean_artifacts
                green_clean{ii}.responsecurve = [];
                green_corrected_clean{ii}.responsecurve = [];
                green_smooth_clean{ii}.responsecurve = [];
            end
            disp('Analyzing green');

        end

         if isfield(fiber, 'isosbestic') 
            iso{ii}.responsecurve = [];
            iso_corrected{ii}.responsecurve = [];
            iso_smooth{ii}.responsecurve = [];
            if fiber.clean_artifacts
                iso_clean{ii}.responsecurve = [];
                iso_corrected_clean{ii}.responsecurve = [];
                iso_smooth_clean{ii}.responsecurve = [];
            end

        end
    
    
        for jj = 1:length(pul)
               
            [~,idx] = min(abs(fiber.timestamps - pul(jj)));
            % Indexes for the window
            idx_range = idx + round(baseline_ints(1)*fiber.sr):idx + round(event_ints(2)*fiber.sr);
            if length(time_vector) +1 == length(idx_range)
                idx_range(end) = [];
            end
            % idx_range = idx + baseline_ints(1)*fiber.sr:idx + event_ints(2)*fiber.sr;
            
            % Verify that window is inside limits. GENERAL COMPUTATIONS
            if min(idx_range) > 0 && max(idx_range) <= length(fiber.timestamps)
                count = count + 1;
                % red
                if isfield(fiber, 'red') 
                    red{ii}.responsecurve = [red{ii}.responsecurve; fiber.red(idx_range)'];
                    red_corrected{jj}.responsecurve = [red_corrected{jj}.responsecurve; fiber.red_corrected(idx_range)'];
                    red_smooth{jj}.responsecurve = [red_smooth{jj}.responsecurve; fiber.red_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        red_clean{jj}.responsecurve = [red_clean{jj}.responsecurve; fiber.red_clean(idx_range)'];
                        red_corrected_clean{jj}.responsecurve = [red_corrected_clean{jj}.responsecurve; fiber.red_corrected_clean(idx_range)'];
                        red_smooth_clean{jj}.responsecurve = [red_smooth_clean{jj}.responsecurve; fiber.red_smooth_clean(idx_range)'];

                    end

                end
    
                % green
                if isfield(fiber, 'green') 
                    green{ii}.responsecurve = [green{ii}.responsecurve; fiber.green(idx_range)'];
                    green_corrected{jj}.responsecurve = [green_corrected{jj}.responsecurve; fiber.green_corrected(idx_range)'];
                    green_smooth{jj}.responsecurve = [green_smooth{jj}.responsecurve; fiber.green_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        green_clean{jj}.responsecurve = [green_clean{jj}.responsecurve; fiber.green_clean(idx_range)'];
                        green_corrected_clean{jj}.responsecurve = [green_corrected_clean{jj}.responsecurve; fiber.green_corrected_clean(idx_range)'];
                        green_smooth_clean{jj}.responsecurve = [green_smooth_clean{jj}.responsecurve; fiber.green_smooth_clean(idx_range)'];
                    end

                end

                % iso
                if isfield(fiber, 'isosbestic') 
                    iso{ii}.responsecurve = [iso{ii}.responsecurve; fiber.isosbestic(idx_range)'];
                    iso_corrected{jj}.responsecurve = [iso_corrected{jj}.responsecurve; fiber.iso_corrected(idx_range)'];
                    iso_smooth{jj}.responsecurve = [iso_smooth{jj}.responsecurve; fiber.iso_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        iso_clean{jj}.responsecurve = [iso_clean{jj}.responsecurve; fiber.iso_clean(idx_range)'];
                        iso_corrected_clean{jj}.responsecurve = [iso_corrected_clean{jj}.responsecurve; fiber.iso_corrected_clean(idx_range)'];
                        iso_smooth_clean{jj}.responsecurve = [iso_smooth_clean{jj}.responsecurve; fiber.iso_smooth_clean(idx_range)'];
                    end

                end
        
                ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{ii}(count) = fiber.timestamps(idx);
            end
        end

        % RED FLUORESCENCE (Ca2+)
    if isfield(fiber, 'red') 
        
        for ii = 1: size(red{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses

                red{jj}.responsecurveZ(ii,:) = (red{jj}.responsecurve(ii,:)...
                    -mean(red{jj}.responsecurve(ii,t_Z)))...
                    ./std(red{jj}.responsecurve(ii,t_Z));

                red_corrected{jj}.responsecurveZ(ii,:) = (red_corrected{jj}.responsecurve(ii,:)...
                    -mean(red_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_corrected{jj}.responsecurve(ii,t_Z));

                red_smooth{jj}.responsecurveZ(ii,:) = (red_smooth{jj}.responsecurve(ii,:)...
                    -mean(red_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_smooth{jj}.responsecurve(ii,t_Z));

                if fiber.clean_artifacts

                    red_clean{jj}.responsecurveZ(ii,:) = (red_clean{jj}.responsecurve(ii,:)...
                    -mean(red_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_clean{jj}.responsecurve(ii,t_Z));

                    red_corrected_clean{jj}.responsecurveZ(ii,:) = (red_corrected_clean{jj}.responsecurve(ii,:)...
                        -mean(red_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(red_corrected_clean{jj}.responsecurve(ii,t_Z));
    
                    red_smooth_clean{jj}.responsecurveZ(ii,:) = (red_smooth_clean{jj}.responsecurve(ii,:)...
                        -mean(red_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(red_smooth_clean{jj}.responsecurve(ii,t_Z));
                   
                end

            end
        end
    end 
    
    % GREEN FLUORESCENCE (eCB)
    if isfield(fiber, 'green') 
        
        for ii = 1:size(green{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses

                green{jj}.responsecurveZ(ii,:) = (green{jj}.responsecurve(ii,:)...
                    -mean(green{jj}.responsecurve(ii,t_Z)))...
                    ./std(green{jj}.responsecurve(ii,t_Z));

                green_corrected{jj}.responsecurveZ(ii,:) = (green_corrected{jj}.responsecurve(ii,:)...
                    -mean(green_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_corrected{jj}.responsecurve(ii,t_Z));

                green_smooth{jj}.responsecurveZ(ii,:) = (green_smooth{jj}.responsecurve(ii,:)...
                    -mean(green_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_smooth{jj}.responsecurve(ii,t_Z));

                if fiber.clean_artifacts
                    green_clean{jj}.responsecurveZ(ii,:) = (green_clean{jj}.responsecurve(ii,:)...
                    -mean(green_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_clean{jj}.responsecurve(ii,t_Z));

                green_corrected_clean{jj}.responsecurveZ(ii,:) = (green_corrected_clean{jj}.responsecurve(ii,:)...
                    -mean(green_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_corrected_clean{jj}.responsecurve(ii,t_Z));

                green_smooth_clean{jj}.responsecurveZ(ii,:) = (green_smooth_clean{jj}.responsecurve(ii,:)...
                    -mean(green_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_smooth_clean{jj}.responsecurve(ii,t_Z));

                end
            end
        end
    end


    % ISOSBESTIC
    if isfield(fiber, 'iso') 
        
        for ii = 1:size(iso{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses
  
                iso{jj}.responsecurveZ(ii,:) = (iso{jj}.responsecurve(ii,:)...
                    -mean(iso{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso{jj}.responsecurve(ii,t_Z));

                iso_corrected{jj}.responsecurveZ(ii,:) = (iso_corrected{jj}.responsecurve(ii,:)...
                    -mean(iso_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_corrected{jj}.responsecurve(ii,t_Z));

                iso_smooth{jj}.responsecurveZ(ii,:) = (iso_smooth{jj}.responsecurve(ii,:)...
                    -mean(iso_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_smooth{jj}.responsecurve(ii,t_Z));

                if fiber.clean_artifacts
                    iso_clean{jj}.responsecurveZ(ii,:) = (iso_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_clean{jj}.responsecurve(ii,t_Z));

                iso_corrected_clean{jj}.responsecurveZ(ii,:) = (iso_corrected_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_corrected_clean{jj}.responsecurve(ii,t_Z));

                iso_smooth_clean{jj}.responsecurveZ(ii,:) = (iso_smooth_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_smooth_clean{jj}.responsecurve(ii,t_Z));

                end
            end
        end
    end
end

else
    % Analyze fiber
    count = 0;
    for ii = 1:nConditions
    
        fprintf('\n Condition %3.i / %3.i \n', ii, nConditions)
    
        pul = pulses.timestamps(pulses.channel == conditions(ii,2) & pulses.duration == conditions(ii,1),1);
        channelPulse = median(pulses.channel(pulses.channel == conditions(ii,2) & pulses.duration == conditions(ii,1)));
        if isempty(pul)
            pul = [0];
        end
        if length(pul)>maxNumberOfPulses
            pul = sort(randsample(pul, maxNumberOfPulses));
        end
    
        if isfield(fiber, 'red') 
            red{ii}.responsecurve = [red{ii}.responsecurve; fiber.red(idx_range)'];
            red_corrected{jj}.responsecurve = [red_corrected{jj}.responsecurve; fiber.red_corrected(idx_range)'];
            red_smooth{jj}.responsecurve = [red_smooth{jj}.responsecurve; fiber.red_smooth(idx_range)'];
            if fiber.clean_artifacts
                red_clean{jj}.responsecurve = [red_clean{jj}.responsecurve; fiber.red_clean(idx_range)'];
                red_corrected_clean{jj}.responsecurve = [red_corrected_clean{jj}.responsecurve; fiber.red_corrected_clean(idx_range)'];
                red_smooth_clean{jj}.responsecurve = [red_smooth_clean{jj}.responsecurve; fiber.red_smooth_clean(idx_range)'];

            end
        end
        if isfield(fiber, 'green') 
            green{ii}.responsecurve = [green{ii}.responsecurve; fiber.green(idx_range)'];
            green_corrected{jj}.responsecurve = [green_corrected{jj}.responsecurve; fiber.green_corrected(idx_range)'];
            green_smooth{jj}.responsecurve = [green_smooth{jj}.responsecurve; fiber.green_smooth(idx_range)'];
            if fiber.clean_artifacts
                green_clean{jj}.responsecurve = [green_clean{jj}.responsecurve; fiber.green_clean(idx_range)'];
                green_corrected_clean{jj}.responsecurve = [green_corrected_clean{jj}.responsecurve; fiber.green_corrected_clean(idx_range)'];
                green_smooth_clean{jj}.responsecurve = [green_smooth_clean{jj}.responsecurve; fiber.green_smooth_clean(idx_range)'];
            end
 
        end

        if isfield(fiber, 'isosbestic') 
            iso{ii}.responsecurve = [iso{ii}.responsecurve; fiber.isosbestic(idx_range)'];
            iso_corrected{jj}.responsecurve = [iso_corrected{jj}.responsecurve; fiber.iso_corrected(idx_range)'];
            iso_smooth{jj}.responsecurve = [iso_smooth{jj}.responsecurve; fiber.iso_smooth(idx_range)'];
            if fiber.clean_artifacts
                iso_clean{jj}.responsecurve = [iso_clean{jj}.responsecurve; fiber.iso_clean(idx_range)'];
                iso_corrected_clean{jj}.responsecurve = [iso_corrected_clean{jj}.responsecurve; fiber.iso_corrected_clean(idx_range)'];
                iso_smooth_clean{jj}.responsecurve = [iso_smooth_clean{jj}.responsecurve; fiber.iso_smooth_clean(idx_range)'];
            end

        end
    
    
        for jj = 1:length(pul)
               
            [~,idx] = min(abs(fiber.timestamps - pul(jj)));
            % Indexes for the window
            idx_range = idx + round(baseline_ints(1)*fiber.sr):idx + round(event_ints(2)*fiber.sr);
            if length(time_vector) +1 == length(idx_range)
                idx_range(end) = [];
            end
            % idx_range = idx + baseline_ints(1)*fiber.sr:idx + event_ints(2)*fiber.sr;
            
            % Verify that window is inside limits. GENERAL COMPUTATIONS
            if min(idx_range) > 0 && max(idx_range) <= length(fiber.timestamps)
                count = count + 1;
                % red
                if isfield(fiber, 'red') 
                    for ii = 1: size(red{jj}.responsecurve,1)
                        if numberOfPulses > minNumberOfPulses
            
                            red{jj}.responsecurveZ(ii,:) = (red{jj}.responsecurve(ii,:)...
                                -mean(red{jj}.responsecurve(ii,t_Z)))...
                                ./std(red{jj}.responsecurve(ii,t_Z));
            
                            red_corrected{jj}.responsecurveZ(ii,:) = (red_corrected{jj}.responsecurve(ii,:)...
                                -mean(red_corrected{jj}.responsecurve(ii,t_Z)))...
                                ./std(red_corrected{jj}.responsecurve(ii,t_Z));
            
                            red_smooth{jj}.responsecurveZ(ii,:) = (red_smooth{jj}.responsecurve(ii,:)...
                                -mean(red_smooth{jj}.responsecurve(ii,t_Z)))...
                                ./std(red_smooth{jj}.responsecurve(ii,t_Z));
            
                            if fiber.clean_artifacts
            
                                red_clean{jj}.responsecurveZ(ii,:) = (red_clean{jj}.responsecurve(ii,:)...
                                -mean(red_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(red_clean{jj}.responsecurve(ii,t_Z));
            
                                red_corrected_clean{jj}.responsecurveZ(ii,:) = (red_corrected_clean{jj}.responsecurve(ii,:)...
                                    -mean(red_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                                    ./std(red_corrected_clean{jj}.responsecurve(ii,t_Z));
                
                                red_smooth_clean{jj}.responsecurveZ(ii,:) = (red_smooth_clean{jj}.responsecurve(ii,:)...
                                    -mean(red_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                                    ./std(red_smooth_clean{jj}.responsecurve(ii,t_Z));
                               
                            end
            
                        end
                    end

                end
    
                % green
                if isfield(fiber, 'green') 
                    for ii = 1:size(green{jj}.responsecurve,1)
                        if numberOfPulses > minNumberOfPulses
            
                            green{jj}.responsecurveZ(ii,:) = (green{jj}.responsecurve(ii,:)...
                                -mean(green{jj}.responsecurve(ii,t_Z)))...
                                ./std(green{jj}.responsecurve(ii,t_Z));
            
                            green_corrected{jj}.responsecurveZ(ii,:) = (green_corrected{jj}.responsecurve(ii,:)...
                                -mean(green_corrected{jj}.responsecurve(ii,t_Z)))...
                                ./std(green_corrected{jj}.responsecurve(ii,t_Z));
            
                            green_smooth{jj}.responsecurveZ(ii,:) = (green_smooth{jj}.responsecurve(ii,:)...
                                -mean(green_smooth{jj}.responsecurve(ii,t_Z)))...
                                ./std(green_smooth{jj}.responsecurve(ii,t_Z));
            
                            if fiber.clean_artifacts
                                green_clean{jj}.responsecurveZ(ii,:) = (green_clean{jj}.responsecurve(ii,:)...
                                -mean(green_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(green_clean{jj}.responsecurve(ii,t_Z));
            
                            green_corrected_clean{jj}.responsecurveZ(ii,:) = (green_corrected_clean{jj}.responsecurve(ii,:)...
                                -mean(green_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(green_corrected_clean{jj}.responsecurve(ii,t_Z));
            
                            green_smooth_clean{jj}.responsecurveZ(ii,:) = (green_smooth_clean{jj}.responsecurve(ii,:)...
                                -mean(green_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(green_smooth_clean{jj}.responsecurve(ii,t_Z));
            
                            end
                        end
                    end

                end

                % iso
                if isfield(fiber, 'isosbestic') 
                    for ii = 1:size(iso{jj}.responsecurve,1)
                        if numberOfPulses > minNumberOfPulses
              
                            iso{jj}.responsecurveZ(ii,:) = (iso{jj}.responsecurve(ii,:)...
                                -mean(iso{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso{jj}.responsecurve(ii,t_Z));
            
                            iso_corrected{jj}.responsecurveZ(ii,:) = (iso_corrected{jj}.responsecurve(ii,:)...
                                -mean(iso_corrected{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso_corrected{jj}.responsecurve(ii,t_Z));
            
                            iso_smooth{jj}.responsecurveZ(ii,:) = (iso_smooth{jj}.responsecurve(ii,:)...
                                -mean(iso_smooth{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso_smooth{jj}.responsecurve(ii,t_Z));
            
                            if fiber.clean_artifacts
                                iso_clean{jj}.responsecurveZ(ii,:) = (iso_clean{jj}.responsecurve(ii,:)...
                                -mean(iso_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso_clean{jj}.responsecurve(ii,t_Z));
            
                            iso_corrected_clean{jj}.responsecurveZ(ii,:) = (iso_corrected_clean{jj}.responsecurve(ii,:)...
                                -mean(iso_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso_corrected_clean{jj}.responsecurve(ii,t_Z));
            
                            iso_smooth_clean{jj}.responsecurveZ(ii,:) = (iso_smooth_clean{jj}.responsecurve(ii,:)...
                                -mean(iso_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                                ./std(iso_smooth_clean{jj}.responsecurve(ii,t_Z));
            
                            end
                        end
                    end

                end


        
                ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{ii}(count) = fiber.timestamps(idx);
            end
        end


        % GREEN FLUORESCENCE (eCB)
    if isfield(fiber, 'green') 
        
        for ii = 1:size(green{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses

                green{jj}.responsecurveZ(ii,:) = (green{jj}.responsecurve(ii,:)...
                    -mean(green{jj}.responsecurve(ii,t_Z)))...
                    ./std(green{jj}.responsecurve(ii,t_Z));

                green_corrected{jj}.responsecurveZ(ii,:) = (green_corrected{jj}.responsecurve(ii,:)...
                    -mean(green_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_corrected{jj}.responsecurve(ii,t_Z));

                green_smooth{jj}.responsecurveZ(ii,:) = (green_smooth{jj}.responsecurve(ii,:)...
                    -mean(green_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_smooth{jj}.responsecurve(ii,t_Z));

                if fiber.clean_artifacts
                    green_clean{jj}.responsecurveZ(ii,:) = (green_clean{jj}.responsecurve(ii,:)...
                    -mean(green_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_clean{jj}.responsecurve(ii,t_Z));

                green_corrected_clean{jj}.responsecurveZ(ii,:) = (green_corrected_clean{jj}.responsecurve(ii,:)...
                    -mean(green_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_corrected_clean{jj}.responsecurve(ii,t_Z));

                green_smooth_clean{jj}.responsecurveZ(ii,:) = (green_smooth_clean{jj}.responsecurve(ii,:)...
                    -mean(green_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_smooth_clean{jj}.responsecurve(ii,t_Z));

                end
            end
        end
    end


    % ISOSBESTIC
    if isfield(fiber, 'iso') 
        
        for ii = 1:size(iso{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses
  
                iso{jj}.responsecurveZ(ii,:) = (iso{jj}.responsecurve(ii,:)...
                    -mean(iso{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso{jj}.responsecurve(ii,t_Z));

                iso_corrected{jj}.responsecurveZ(ii,:) = (iso_corrected{jj}.responsecurve(ii,:)...
                    -mean(iso_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_corrected{jj}.responsecurve(ii,t_Z));

                iso_smooth{jj}.responsecurveZ(ii,:) = (iso_smooth{jj}.responsecurve(ii,:)...
                    -mean(iso_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_smooth{jj}.responsecurve(ii,t_Z));

                if fiber.clean_artifacts
                    iso_clean{jj}.responsecurveZ(ii,:) = (iso_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_clean{jj}.responsecurve(ii,t_Z));

                iso_corrected_clean{jj}.responsecurveZ(ii,:) = (iso_corrected_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_corrected_clean{jj}.responsecurve(ii,t_Z));

                iso_smooth_clean{jj}.responsecurveZ(ii,:) = (iso_smooth_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_smooth_clean{jj}.responsecurve(ii,t_Z));

                end
            end
        end
    end
end

end

  
if plt

    if length(blocks) > 1

        % Green
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(green{ii}.responsecurveZ), green{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
    
            zmean = mean(green{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

        end

        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end

        cmap = jet(length(blocks));


        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, green{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end

        % Green corrected
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(green_corrected{ii}.responsecurveZ), green_corrected{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
    
            zmean = mean(green_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

        end

        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end

        cmap = jet(length(blocks));


        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, green_corrected{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end


        % Green smooth
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(green_smooth{ii}.responsecurveZ), green_smooth{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
    
            zmean = mean(green_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

        end

        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end

        cmap = jet(length(blocks));


        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, green_smooth{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end

        if fiber.clean_artifacts
            % Green clean
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_clean{ii}.responsecurveZ), green_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            end
    
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
    
            cmap = jet(length(blocks));
    
    
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
    
            % Green corrected
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_corrected_clean{ii}.responsecurveZ), green_corrected_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_corrected_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            end
    
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
    
            cmap = jet(length(blocks));
    
    
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_corrected_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
    
    
            % Green smooth
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_smooth_clean{ii}.responsecurveZ), green_smooth_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_smooth_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            end
    
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_smooth_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
    
            cmap = jet(length(blocks));
    
    
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_smooth_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_smooth_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end

        end

        % red
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(red{ii}.responsecurveZ), red{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(red{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, red{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        % red corrected
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(red_corrected{ii}.responsecurveZ), red_corrected{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(red_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, red_corrected{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        
        % red smooth
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(red_smooth{ii}.responsecurveZ), red_smooth{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(red_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, red_smooth{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_red_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        if fiber.clean_artifacts
            % red clean
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(red_clean{ii}.responsecurveZ), red_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(red_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(red_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, red_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
            % red corrected
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(red_corrected_clean{ii}.responsecurveZ), red_corrected_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(red_corrected_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(red_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, red_corrected_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
        
            % red smooth
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(red_smooth_clean{ii}.responsecurveZ), red_smooth_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(red_smooth_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(red_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_smooth_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, red_smooth_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_red_smooth_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
        end


        % iso
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(iso{ii}.responsecurveZ), iso{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(iso{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, iso{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        % iso corrected
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(iso_corrected{ii}.responsecurveZ), iso_corrected{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(iso_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, iso_corrected{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_corrected_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        
        % iso smooth
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for ii = 1:length(blocks)    
            subplot(2,ceil(size(blocks,1)/2),ii);  
                
            imagesc(time_vector, 1:size(iso_smooth{ii}.responsecurveZ), iso_smooth{ii}.responsecurveZ); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            caxis([-c_axis c_axis]);
            set(gca,'YDir','normal');
            hold on; 
        
            zmean = mean(iso_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        end
        
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
        end
        
        cmap = jet(length(blocks));
        
        
        figure;
        for ii = 1:length(blocks)
            plotFill(time_vector, iso_smooth{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
        end
        if savePlot
            saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
        end
        
        if fiber.clean_artifacts
            % iso clean
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(iso_clean{ii}.responsecurveZ), iso_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(iso_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(iso_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, iso_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
            % iso corrected
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(iso_corrected_clean{ii}.responsecurveZ), iso_corrected_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(iso_corrected_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(iso_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, iso_corrected_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_corrected_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
        
            % iso smooth
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(iso_smooth_clean{ii}.responsecurveZ), iso_smooth_clean{ii}.responsecurveZ); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(iso_smooth_clean{ii}.responsecurveZ);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(iso_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            end
        
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_smooth_clean_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end
        
            cmap = jet(length(blocks));
        
        
            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, iso_smooth_clean{ii}.responsecurveZ,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_smooth_clean',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end
        
        end

    else

            for ii = 1:nConditions
    
               
        
        
            end
    end

end



fiberOptogeneticResponse = [];


    


 





end