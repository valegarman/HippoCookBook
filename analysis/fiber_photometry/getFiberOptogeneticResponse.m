function [fiberOptogeneticResponse] = getFiberOptogeneticResponse(varargin)

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
    fiber = getSessionFiberPhotometry_temp('force',reload_fiber);
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
t_Z = time_vector > -3 & time_vector < -1;

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
            red_smooth{ii}.responsecurve = [];
            red_normalized{ii}.responsecurve = [];
            red_PP{ii}.responsecurve = [];
            % disp('Processing red channel')
        end
        if isfield(fiber, 'green') 
            green{ii}.responsecurve = [];
            green_smooth{ii}.responsecurve = [];
            green_normalized{ii}.responsecurve = [];
            green_PP{ii}.responsecurve = [];
            % disp('Processing green channel')
        end

         if isfield(fiber, 'isosbestic') 
            iso{ii}.responsecurve = [];
            iso_smooth{ii}.responsecurve = [];
            iso_normalized{ii}.responsecurve = [];
            iso_PP{ii}.responsecurve = [];
            % disp('Processing isosbestic channel')
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
                    red_smooth{ii}.responsecurve = [red_smooth{ii}.responsecurve; fiber.red_fpa.fSmoothed(idx_range)'];
                    red_normalized{ii}.responsecurve = [red_normalized{ii}.responsecurve; fiber.red_fpa.fNormalized(idx_range)'];
    
                    red_PP{ii}.responsecurve = [red_PP{ii}.responsecurve; fiber.red_PP.red_dFF_Smoothed(idx_range)'];
                end
    
                % green
                if isfield(fiber, 'green') 
                    green{ii}.responsecurve = [green{ii}.responsecurve; fiber.green(idx_range)'];
                    green_smooth{ii}.responsecurve = [green_smooth{ii}.responsecurve; fiber.green_fpa.fSmoothed(idx_range)'];
                    green_normalized{ii}.responsecurve = [green_normalized{ii}.responsecurve; fiber.green_fpa.fNormalized(idx_range)'];
    
                    green_PP{ii}.responsecurve = [green_PP{ii}.responsecurve; fiber.green_PP.green_dFF_Smoothed(idx_range)'];
                end

                % iso
                if isfield(fiber, 'isosbestic') 
                    iso{ii}.responsecurve = [iso{ii}.responsecurve; fiber.isosbestic(idx_range)'];
                    % iso_smooth{ii}.responsecurve = [iso_smooth{ii}.responsecurve; fiber.isosbestic_fpa.fSmoothed(idx_range)'];
                    % iso_normalized{ii}.responsecurve = [iso_normalized{ii}.responsecurve; fiber.isosbestic_fpa.fNormalized(idx_range)'];
    
                    iso_PP{ii}.responsecurve = [iso_PP{ii}.responsecurve; fiber.iso_PP.iso_dFF_Smoothed(idx_range)'];
                end
        
                % ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{ii}(count) = fiber.timestamps(idx);
            end
        end
    
    
        % t_duringPulse = time_vector > event_ints(1) & time_vector < event_ints(2);
        % t_beforePulse = time_vector > baseline_ints(1) & time_vector < baseline_ints(2);
        % t_Z = time_vector <= event_ints(1);
    
        % GREEN FLUORESCENCE (eCB)
        if isfield(fiber, 'green') 
            % green
            f_green_normalized_prctl20 = prctile(fiber.green_fpa.fNormalized,20);
            
            for jj = 1:size(green{ii}.responsecurve,1)
                green{ii}.responsecurveSmooth(jj,:) = smooth(green{ii}.responsecurve(jj,:));
                green{ii}.responsecurveZ(jj,:) = (green{ii}.responsecurve(jj,:)...
                    -mean(green{ii}.responsecurve(jj,t_Z)))...
                    /std(green{ii}.responsecurve(jj,t_Z));
                green{ii}.responsecurveZSmooth(jj,:) = smooth(green{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green{ii}.prctile(jj,:) = green{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
        
            end
            
            
            % green smooth
            f_green_smooth_prctl20 = prctile(fiber.green_fpa.fSmoothed,20);
            
            for jj = 1:size(green_smooth{ii}.responsecurve,1)
                green_smooth{ii}.responsecurveSmooth(jj,:) = smooth(green_smooth{ii}.responsecurve(jj,:));
                green_smooth{ii}.responsecurveZ(jj,:) = (green_smooth{ii}.responsecurve(jj,:)...
                    -mean(green_smooth{ii}.responsecurve(jj,t_Z)))...
                    /std(green_smooth{ii}.responsecurve(jj,t_Z));
                green_smooth{ii}.responsecurveZSmooth(jj,:) = smooth(green_smooth{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_smooth{ii}.prctile(jj,:) = green_smooth{ii}.responsecurve(jj,:) - f_green_smooth_prctl20/f_green_smooth_prctl20;
            end
            
            
            % green normalized
            f_green_normalized_prctl20 = prctile(fiber.green_fpa.fNormalized,20);
            
            for jj = 1:size(green_normalized{ii}.responsecurve,1)
                green_normalized{ii}.responsecurveSmooth(jj,:) = smooth(green_normalized{ii}.responsecurve(jj,:));
                green_normalized{ii}.responsecurveZ(jj,:) = (green_normalized{ii}.responsecurve(jj,:)...
                    -mean(green_normalized{ii}.responsecurve(jj,t_Z)))...
                    /std(green_normalized{ii}.responsecurve(jj,t_Z));
                green_normalized{ii}.responsecurveZSmooth(jj,:) = smooth(green_normalized{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_normalized{ii}.prctile(jj,:) = green_normalized{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
            end
        
            % green PP
            f_green_PP_prctl20 = prctile(fiber.green_PP.green_dFF_Smoothed,20);
            
            for jj = 1:size(green_PP{ii}.responsecurve,1)
                green_PP{ii}.responsecurveSmooth(jj,:) = smooth(green_PP{ii}.responsecurve(jj,:));
                green_PP{ii}.responsecurveZ(jj,:) = (green_PP{ii}.responsecurve(jj,:)...
                    -mean(green_PP{ii}.responsecurve(jj,t_Z)))...
                    /std(green_PP{ii}.responsecurve(jj,t_Z));
                green_PP{ii}.responsecurveZSmooth(jj,:) = smooth(green_PP{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_PP{ii}.prctile(jj,:) = green_PP{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
            end
        end

        % ISOSBESTIC (CONTROL)
        if isfield(fiber, 'isosbestic') 
            % ISO
            % f_iso_normalized_prctl20 = prctile(fiber.iso.fNormalized,20);
            
            for jj = 1:size(iso{ii}.responsecurve,1)
                iso{ii}.responsecurveSmooth(jj,:) = smooth(iso{ii}.responsecurve(jj,:));
                iso{ii}.responsecurveZ(jj,:) = (iso{ii}.responsecurve(jj,:)...
                    -mean(iso{ii}.responsecurve(jj,t_Z)))...
                    /std(iso{ii}.responsecurve(jj,t_Z));
                iso{ii}.responsecurveZSmooth(jj,:) = smooth(iso{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                % iso{ii}.prctile(jj,:) = iso{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            
        
            end
            
            
            % ISO smooth
            % f_iso_smooth_prctl20 = prctile(fiber.iso_fpa.fSmoothed,20);
            % 
            % for jj = 1:size(iso_smooth{ii}.responsecurve,1)
            %     iso_smooth{ii}.responsecurveSmooth(jj,:) = smooth(iso_smooth{ii}.responsecurve(jj,:));
            %     iso_smooth{ii}.responsecurveZ(jj,:) = (iso_smooth{ii}.responsecurve(jj,:)...
            %         -mean(iso_smooth{ii}.responsecurve(jj,t_Z)))...
            %         /std(iso_smooth{ii}.responsecurve(jj,t_Z));
            %     iso_smooth{ii}.responsecurveZSmooth(jj,:) = smooth(iso_smooth{ii}.responsecurveZ(jj,:));
            % 
            %     % Using the 20th percentile
            %     iso_smooth{ii}.prctile(jj,:) = iso_smooth{ii}.responsecurve(jj,:) - f_iso_smooth_prctl20/f_iso_smooth_prctl20;
            % end
            
            
            % ISO normalized
            % f_iso_normalized_prctl20 = prctile(fiber.iso_fpa.fNormalized,20);
            % 
            % for jj = 1:size(iso_normalized{ii}.responsecurve,1)
            %     iso_normalized{ii}.responsecurveSmooth(jj,:) = smooth(iso_normalized{ii}.responsecurve(jj,:));
            %     iso_normalized{ii}.responsecurveZ(jj,:) = (iso_normalized{ii}.responsecurve(jj,:)...
            %         -mean(iso_normalized{ii}.responsecurve(jj,t_Z)))...
            %         /std(iso_normalized{ii}.responsecurve(jj,t_Z));
            %     iso_normalized{ii}.responsecurveZSmooth(jj,:) = smooth(iso_normalized{ii}.responsecurveZ(jj,:));
            % 
            %     % Using the 20th percentile
            %     iso_normalized{ii}.prctile(jj,:) = iso_normalized{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            % 
            % end
        
            % ISO PP
            % f_iso_PP_prctl20 = prctile(fiber.iso_PP.iso_dFF_Smoothed,20);
            
            for jj = 1:size(iso_PP{ii}.responsecurve,1)
                iso_PP{ii}.responsecurveSmooth(jj,:) = smooth(iso_PP{ii}.responsecurve(jj,:));
                iso_PP{ii}.responsecurveZ(jj,:) = (iso_PP{ii}.responsecurve(jj,:)...
                    -mean(iso_PP{ii}.responsecurve(jj,t_Z)))...
                    /std(iso_PP{ii}.responsecurve(jj,t_Z));
                iso_PP{ii}.responsecurveZSmooth(jj,:) = smooth(iso_PP{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                % iso_PP{ii}.prctile(jj,:) = iso_PP{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            
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
            red{ii}.responsecurve = [];
            red_smooth{ii}.responsecurve = [];
            red_normalized{ii}.responsecurve = [];
            red_PP{ii}.responsecurve = [];
            % disp('Processing red channel')
        end
        if isfield(fiber, 'green') 
            green{ii}.responsecurve = [];
            green_smooth{ii}.responsecurve = [];
            green_normalized{ii}.responsecurve = [];
            green_PP{ii}.responsecurve = [];
            % disp('Processing green channel')
        end

        if isfield(fiber, 'isosbestic') 
            iso{ii}.responsecurve = [];
            iso_smooth{ii}.responsecurve = [];
            iso_normalized{ii}.responsecurve = [];
            iso_PP{ii}.responsecurve = [];
            % disp('Processing isosbestic channel')
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
                    red_smooth{ii}.responsecurve = [red_smooth{ii}.responsecurve; fiber.red_fpa.fSmoothed(idx_range)'];
                    red_normalized{ii}.responsecurve = [red_normalized{ii}.responsecurve; fiber.red_fpa.fNormalized(idx_range)'];
    
                    red_PP{ii}.responsecurve = [red_PP{ii}.responsecurve; fiber.red_PP.red_dFF_Smoothed(idx_range)'];
                end
    
                % green
                if isfield(fiber, 'green') 
                    green{ii}.responsecurve = [green{ii}.responsecurve; fiber.green(idx_range)'];
                    green_smooth{ii}.responsecurve = [green_smooth{ii}.responsecurve; fiber.green_fpa.fSmoothed(idx_range)'];
                    green_normalized{ii}.responsecurve = [green_normalized{ii}.responsecurve; fiber.green_fpa.fNormalized(idx_range)'];
    
                    green_PP{ii}.responsecurve = [green_PP{ii}.responsecurve; fiber.green_PP.green_dFF_Smoothed(idx_range)'];
                end

                % iso
                if isfield(fiber, 'isosbestic') 
                    iso{ii}.responsecurve = [iso{ii}.responsecurve; fiber.isosbestic(idx_range)'];
                    % iso_smooth{ii}.responsecurve = [iso_smooth{ii}.responsecurve; fiber.isosbestic_fpa.fSmoothed(idx_range)'];
                    % iso_normalized{ii}.responsecurve = [iso_normalized{ii}.responsecurve; fiber.isosbestic.fNormalized(idx_range)'];
    
                    iso_PP{ii}.responsecurve = [iso_PP{ii}.responsecurve; fiber.iso_PP.iso_dFF_Smoothed(idx_range)'];
                end


        
                % ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{ii}(count) = fiber.timestamps(idx);
            end
        end
    
    
        t_duringPulse = time_vector > event_ints(1) & time_vector < event_ints(2);
        t_beforePulse = time_vector > baseline_ints(1) & time_vector < baseline_ints(2);
        t_Z = time_vector <= event_ints(1);
    
        % GREEN FLUORESCENCE (eCB)
        if isfield(fiber, 'green') 
            % green
            f_green_normalized_prctl20 = prctile(fiber.green_fpa.fNormalized,20);
            
            for jj = 1:size(green{ii}.responsecurve,1)
                green{ii}.responsecurveSmooth(jj,:) = smooth(green{ii}.responsecurve(jj,:));
                green{ii}.responsecurveZ(jj,:) = (green{ii}.responsecurve(jj,:)...
                    -mean(green{ii}.responsecurve(jj,t_Z)))...
                    /std(green{ii}.responsecurve(jj,t_Z));
                green{ii}.responsecurveZSmooth(jj,:) = smooth(green{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green{ii}.prctile(jj,:) = green{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
        
            end
            
            
            % green smooth
            f_green_smooth_prctl20 = prctile(fiber.green_fpa.fSmoothed,20);
            
            for jj = 1:size(green_smooth{ii}.responsecurve,1)
                green_smooth{ii}.responsecurveSmooth(jj,:) = smooth(green_smooth{ii}.responsecurve(jj,:));
                green_smooth{ii}.responsecurveZ(jj,:) = (green_smooth{ii}.responsecurve(jj,:)...
                    -mean(green_smooth{ii}.responsecurve(jj,t_Z)))...
                    /std(green_smooth{ii}.responsecurve(jj,t_Z));
                green_smooth{ii}.responsecurveZSmooth(jj,:) = smooth(green_smooth{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_smooth{ii}.prctile(jj,:) = green_smooth{ii}.responsecurve(jj,:) - f_green_smooth_prctl20/f_green_smooth_prctl20;
            end
            
            
            % green normalized
            f_green_normalized_prctl20 = prctile(fiber.green_fpa.fNormalized,20);
            
            for jj = 1:size(green_normalized{ii}.responsecurve,1)
                green_normalized{ii}.responsecurveSmooth(jj,:) = smooth(green_normalized{ii}.responsecurve(jj,:));
                green_normalized{ii}.responsecurveZ(jj,:) = (green_normalized{ii}.responsecurve(jj,:)...
                    -mean(green_normalized{ii}.responsecurve(jj,t_Z)))...
                    /std(green_normalized{ii}.responsecurve(jj,t_Z));
                green_normalized{ii}.responsecurveZSmooth(jj,:) = smooth(green_normalized{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_normalized{ii}.prctile(jj,:) = green_normalized{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
            end
        
            % green PP
            f_green_PP_prctl20 = prctile(fiber.green_PP.green_dFF_Smoothed,20);
            
            for jj = 1:size(green_PP{ii}.responsecurve,1)
                green_PP{ii}.responsecurveSmooth(jj,:) = smooth(green_PP{ii}.responsecurve(jj,:));
                green_PP{ii}.responsecurveZ(jj,:) = (green_PP{ii}.responsecurve(jj,:)...
                    -mean(green_PP{ii}.responsecurve(jj,t_Z)))...
                    /std(green_PP{ii}.responsecurve(jj,t_Z));
                green_PP{ii}.responsecurveZSmooth(jj,:) = smooth(green_PP{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                green_PP{ii}.prctile(jj,:) = green_PP{ii}.responsecurve(jj,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
            
            end
        end


        % ISOSBESTIC (CONTROL)
        if isfield(fiber, 'isosbestic') 
            % ISO
            f_iso_normalized_prctl20 = prctile(fiber.iso.fNormalized,20);
            
            for jj = 1:size(iso{ii}.responsecurve,1)
                iso{ii}.responsecurveSmooth(jj,:) = smooth(iso{ii}.responsecurve(jj,:));
                iso{ii}.responsecurveZ(jj,:) = (iso{ii}.responsecurve(jj,:)...
                    -mean(iso{ii}.responsecurve(jj,t_Z)))...
                    /std(iso{ii}.responsecurve(jj,t_Z));
                iso{ii}.responsecurveZSmooth(jj,:) = smooth(iso{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                iso{ii}.prctile(jj,:) = iso{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            
        
            end
            
            
            % ISO smooth
            % f_iso_smooth_prctl20 = prctile(fiber.iso_fpa.fSmoothed,20);
            % 
            % for jj = 1:size(iso_smooth{ii}.responsecurve,1)
            %     iso_smooth{ii}.responsecurveSmooth(jj,:) = smooth(iso_smooth{ii}.responsecurve(jj,:));
            %     iso_smooth{ii}.responsecurveZ(jj,:) = (iso_smooth{ii}.responsecurve(jj,:)...
            %         -mean(iso_smooth{ii}.responsecurve(jj,t_Z)))...
            %         /std(iso_smooth{ii}.responsecurve(jj,t_Z));
            %     iso_smooth{ii}.responsecurveZSmooth(jj,:) = smooth(iso_smooth{ii}.responsecurveZ(jj,:));
            % 
            %     % Using the 20th percentile
            %     iso_smooth{ii}.prctile(jj,:) = iso_smooth{ii}.responsecurve(jj,:) - f_iso_smooth_prctl20/f_iso_smooth_prctl20;
            % end
            
            
            % ISO normalized
            % f_iso_normalized_prctl20 = prctile(fiber.iso_fpa.fNormalized,20);
            % 
            % for jj = 1:size(iso_normalized{ii}.responsecurve,1)
            %     iso_normalized{ii}.responsecurveSmooth(jj,:) = smooth(iso_normalized{ii}.responsecurve(jj,:));
            %     iso_normalized{ii}.responsecurveZ(jj,:) = (iso_normalized{ii}.responsecurve(jj,:)...
            %         -mean(iso_normalized{ii}.responsecurve(jj,t_Z)))...
            %         /std(iso_normalized{ii}.responsecurve(jj,t_Z));
            %     iso_normalized{ii}.responsecurveZSmooth(jj,:) = smooth(iso_normalized{ii}.responsecurveZ(jj,:));
            % 
            %     % Using the 20th percentile
            %     iso_normalized{ii}.prctile(jj,:) = iso_normalized{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            % 
            % end
        
            % ISO PP
            f_iso_PP_prctl20 = prctile(fiber.iso_PP.iso_dFF_Smoothed,20);
            
            for jj = 1:size(iso_PP{ii}.responsecurve,1)
                iso_PP{ii}.responsecurveSmooth(jj,:) = smooth(iso_PP{ii}.responsecurve(jj,:));
                iso_PP{ii}.responsecurveZ(jj,:) = (iso_PP{ii}.responsecurve(jj,:)...
                    -mean(iso_PP{ii}.responsecurve(jj,t_Z)))...
                    /std(iso_PP{ii}.responsecurve(jj,t_Z));
                iso_PP{ii}.responsecurveZSmooth(jj,:) = smooth(iso_PP{ii}.responsecurveZ(jj,:));
        
                % Using the 20th percentile
                iso_PP{ii}.prctile(jj,:) = iso_PP{ii}.responsecurve(jj,:) - f_iso_normalized_prctl20/f_iso_normalized_prctl20;
            
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
                    
                imagesc(time_vector, 1:size(green{ii}.responsecurveZSmooth), green{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end


            % Green normalized
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_normalized{ii}.responsecurveZSmooth), green_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_normalized{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_normalized{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_normalized_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_normalized{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_normalized_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end


            % Green smooth
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_smooth{ii}.responsecurveZSmooth), green_smooth{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_smooth{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_smooth{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_smooth{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end



            % Green_PP
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(green_PP{ii}.responsecurveZSmooth), green_PP{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(green_PP{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(green_PP{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_PP_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, green_PP{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_green_PP_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end



            % ISO
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(iso{ii}.responsecurveZSmooth), iso{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(iso{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(iso{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, iso{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end


            % ISO normalized
            % figure;
            % set(gcf,'Position',[100 -100 2500 1200])
            % for ii = 1:length(blocks)    
            %     subplot(2,ceil(size(blocks,1)/2),ii);  
            % 
            %     imagesc(time_vector, 1:size(iso_normalized{ii}.responsecurveZSmooth), iso_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
            %     colormap(jet);  % Código de colores para visualizar cambios en la actividad
            %     colorbar;  % Agregar barra de color
            %     xlabel('Time (s)');
            %     ylabel(['Trials ', eventType]);
            %     title(['eCB during ', eventType]);
            %     caxis([-c_axis c_axis]);
            %     set(gca,'YDir','normal');
            %     hold on; 
            % 
            %     zmean = mean(iso_normalized{ii}.responsecurveZSmooth);
            %     zmean = zmean-min(zmean); 
            %     zmean = zmean/max(zmean) * size(iso_normalized{ii}.responsecurveZSmooth,1)+1 * std(zmean);
            %     plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            %     xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
            % 
            % end
            % 
            % if savePlot
            %     saveas(gca,['SummaryFigures\fiber_iso_normalized_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            % end
            % 
            % cmap = jet(length(blocks));
            % 
            % 
            % figure;
            % for ii = 1:length(blocks)
            %     plotFill(time_vector, iso_normalized{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            % end
            % if savePlot
            %     saveas(gca,['SummaryFigures\fiber_iso_normalized_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            % end


            % ISO smooth
            % figure;
            % set(gcf,'Position',[100 -100 2500 1200])
            % for ii = 1:length(blocks)    
            %     subplot(2,ceil(size(blocks,1)/2),ii);  
            % 
            %     imagesc(time_vector, 1:size(iso_smooth{ii}.responsecurveZSmooth), iso_smooth{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
            %     colormap(jet);  % Código de colores para visualizar cambios en la actividad
            %     colorbar;  % Agregar barra de color
            %     xlabel('Time (s)');
            %     ylabel(['Trials ', eventType]);
            %     title(['eCB during ', eventType]);
            %     caxis([-c_axis c_axis]);
            %     set(gca,'YDir','normal');
            %     hold on; 
            % 
            %     zmean = mean(iso_smooth{ii}.responsecurveZSmooth);
            %     zmean = zmean-min(zmean); 
            %     zmean = zmean/max(zmean) * size(iso_smooth{ii}.responsecurveZSmooth,1)+1 * std(zmean);
            %     plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            %     xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
            % 
            % end
            % 
            % if savePlot
            %     saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            % end
            % 
            % cmap = jet(length(blocks));
            % 
            % 
            % figure;
            % for ii = 1:length(blocks)
            %     plotFill(time_vector, iso_smooth{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            % end
            % if savePlot
            %     saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            % end



            % ISO_PP
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            for ii = 1:length(blocks)    
                subplot(2,ceil(size(blocks,1)/2),ii);  
                    
                imagesc(time_vector, 1:size(iso_PP{ii}.responsecurveZSmooth), iso_PP{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                colormap(jet);  % Código de colores para visualizar cambios en la actividad
                colorbar;  % Agregar barra de color
                xlabel('Time (s)');
                ylabel(['Trials ', eventType]);
                title(['eCB during ', eventType]);
                caxis([-c_axis c_axis]);
                set(gca,'YDir','normal');
                hold on; 
        
                zmean = mean(iso_PP{ii}.responsecurveZSmooth);
                zmean = zmean-min(zmean); 
                zmean = zmean/max(zmean) * size(iso_PP{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

            end

            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_PP_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks.png']);
            end

            cmap = jet(length(blocks));


            figure;
            for ii = 1:length(blocks)
                plotFill(time_vector, iso_PP{ii}.responsecurveZSmooth,'color',cmap(ii,:),'smoothOpt',2);
            end
            if savePlot
                saveas(gca,['SummaryFigures\fiber_iso_PP_',eventType,'_',num2str(length(blocks)),'_stimulation_blocks_mean.png']);
            end

            % eCB vs ISO
            figure;
            plotFill(time_vector,iso{ii}.responsecurveZSmooth,'color',[163 73 164]/255);
            hold on;
            plotFill(time_vector,green{ii}.responsecurveZSmooth,'color',[0 1 0]);
            % ylim([-1 1]);
            xlim([-3 3])
            saveas(gca,['SummaryFigures\fiber_iso_vs_green_',eventType,'.png']);
    
            figure;
            plotFill(time_vector,iso_PP{ii}.responsecurveZSmooth,'color',[163 73 164]/255);
            hold on;
            plotFill(time_vector,green_PP{ii}.responsecurveZSmooth,'color',[0 1 0]);
            % ylim([-1 1]);
            xlim([-3 3])
            saveas(gca,['SummaryFigures\fiber_iso_vs_green_PP_',eventType,'.png']);




        else
    
                for ii = 1:nConditions
            
                    
            
        
                    if isfield(fiber,'green')
        
        
                        figure;
                        imagesc(time_vector, 1:size(green{ii}.responsecurveZSmooth,1), green{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                        colormap(jet);  % Código de colores para visualizar cambios en la actividad
                        colorbar;  % Agregar barra de color
        
        
                        figure;
                        imagesc(time_vector, 1:size(green{ii}.responsecurveZSmooth), green{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                        colormap(jet);  % Código de colores para visualizar cambios en la actividad
                        colorbar;  % Agregar barra de color
                        xlabel('Time (s)');
                        ylabel(['Trials ', eventType]);
                        title(['eCB during ', eventType]);
                        caxis([-c_axis c_axis]);
                        set(gca,'YDir','normal');
                        hold on; 
                
                        zmean = mean(green{ii}.responsecurveZSmooth);
                        zmean = zmean-min(zmean); 
                        zmean = zmean/max(zmean) * size(green{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                        plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                        xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
        
        
        
                        figure;
                        imagesc(time_vector, 1:size(green{ii}.responsecurveZSmooth,1), green_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                        colormap(jet);  % Código de colores para visualizar cambios en la actividad
                        colorbar;  % Agregar barra de color
        
        
                        figure;
                        imagesc(time_vector, 1:size(green_normalized{ii}.responsecurveZSmooth), green_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                        colormap(jet);  % Código de colores para visualizar cambios en la actividad
                        colorbar;  % Agregar barra de color
                        xlabel('Time (s)');
                        ylabel(['Trials ', eventType]);
                        title(['eCB during ', eventType]);
                        caxis([-c_axis c_axis]);
                        set(gca,'YDir','normal');
                        hold on; 
                
                        zmean = mean(green_normalized{ii}.responsecurveZSmooth);
                        zmean = zmean-min(zmean); 
                        zmean = zmean/max(zmean) * size(green_normalized{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                        plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                        xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
                
                        if savePlot
                            if restrict_fiber_epochs | ~isempty(savePlotAs)
                                saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',save_plt_as{ii},'.png']);
                            else
                                saveas(gca,['SummaryFigures\fiber_green_',eventType,'.png']);
                            end
                        end
                
                        
                        figure;
                        imagesc(time_vector, 1:count, green_PP{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
                        colormap(jet);  % Código de colores para visualizar cambios en la actividad
                        colorbar;  % Agregar barra de color
                        xlabel('Time (s)');
                        ylabel(['Trials ', eventType]);
                        title(['eCB during ', eventType]);
                        caxis([-c_axis c_axis]);
                        set(gca,'YDir','normal');
                        hold on; 
                
                        zmean = mean(green_PP{ii}.responsecurveZSmooth);
                        zmean = zmean-min(zmean); 
                        zmean = zmean/max(zmean) * size(green_PP{ii}.responsecurveZSmooth,1)+1 * std(zmean);
                        plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
                        xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
                
                        if savePlot
                            if restrict_fiber_epochs | ~isempty(savePlotAs)
                                saveas(gca,['SummaryFigures\fiber_green_PP_',eventType,'_',save_plt_as{ii},'.png']);
                            else
                                saveas(gca,['SummaryFigures\fiber_green_PP_',eventType,'.png']);
                            end
                        end
            
            
                        figure;
                        plotFill(time_vector,green_PP{ii}.responsecurveZSmooth,'color',[0 1 0]);
                        ylim([-1 1]);
            
                        if savePlot
                            if restrict_fiber_epochs | ~isempty(savePlotAs)
                                saveas(gca,['SummaryFigures\fiber_green_mean_PP_',eventType,'_',save_plt_as{ii},'.png']);
                            else
                                saveas(gca,['SummaryFigures\fiber_green_mean_PP_',eventType,'.png']);
                            end
                        end
            
            
            
                    end
            
                    % if savePlot
                    %     if restrict_fiber_epochs | ~isempty(savePlotAs)
                    %         saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',save_plt_as{ii},'.png']);
                    %     else
                    %         saveas(gca,['SummaryFigures\fiber_green_',eventType,'.png']);
                    %     end
                    % end
            
            
            
                end
        end

    end



fiberOptogeneticResponse = [];


    


 





end