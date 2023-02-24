function [spkEventTimes] = getSpikesRank(varargin)
% [SpkEventTimes] = bz_getSpikesRank()
% Saves spike times of different units in different ways:
%   1. Absolute and relative time of spikes by unit and by event
%   2. Absolute and relative time of spikes by unit
%   3. Absolute and relative time of spikes by event
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    full path where session is located (default pwd)
%                   e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%     'events'      It can be either the following options:
%                   1. Buzcode structure for specific events (ripples, UDStates, ...)
%                      By default it will load ripples (output from bz_DetectSWR).
%                      Specifically, if not provided, it loads this event 
%                      structure from 'basepath' (if provided), or from current
%                      folder (if not). Its internal structure must have the 
%                      following field:
%                        .timestamps: Nx2 matrix with starting and ending times 
%                                     (in segs) of each event.
%                   2. A Nx2 matrix with starting and ending times (in segs)
%                      of each event, just like the .timestamps field.
%                   (N: number of events)
%                   3. Event type to compute rank from:
%                        - (upstates) UP states (DOWN TO DOWN sequences),
%                        from UDStates.events.mat
%                        - (ripples) Ripples from ripples.events.mat
%     'spikes'      buzcode event structure (from bz_GetSpikes). 
%                   If not provided, it loads it from 'basepath' (if provided),
%                   or from current folder (if not)
%     'UIDs'        A Mx1 boolean matrix with 1s for units to be considered
%                   and 0s for units to be discarded.
%                   (M: number of units)
%     'padding'     extra time after event end to still search for spikes. 
%                   (default is 0.05 seg)
%     'saveMat'   	Saves file, logical (default: true) 
%
%    =========================================================================
%
% OUTPUTS
%
% spkEventTimes structure with the followin fields:
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%	  .UnitEventAbs	 MxN cell matrix. In each cell, absolute times of
%                    spikes for that particular unit and event
%	  .UnitEventRel	 MxN cell matrix. In each cell, relative times of
%                    spikes (relative to the starting time of events) for 
%                    that particular unit and event
%	  .UnitAbs		 1xM cell matrix. In each cell, absolute times of
%                    spikes for that particular unit across all events
%	  .UnitRel		 1xM cell matrix. In each cell, relative times of
%                    spikes for that particular unit across all events
%	  .EventAbs		 3xN cell matrix. In the first row, absolute times of
%                    spikes for that particular event across all units. In
%                    the second row, the UID associated to the above spike.
%                    In third row, the position within the UID vector of
%                    the above spike.
%	  .EventRel		 3xN cell matrix. In the first row, relative times of
%                    spikes for that particular event across all units. In
%                    the second row, the UID associated to the above spike.
%                    In third row, the position within the UID vector of
%                    the above spike.
%
%    =========================================================================
%
%  See also bz_RankOrder
%
%
%
%    Antonio FR, 2017
% Convert to buzcode format: Andrea Navas-Olive, 2019

% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'events',[], @(x) isnumeric(x) || isstruct(x) || ischar(x));
addParameter(p,'spikes',{},@isstruct);
addParameter(p,'UIDs',[],@islogical);
addParameter(p,'padding',0.05,@isnumeric);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'includeIntervals',[],@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
events = p.Results.events;
spikes = p.Results.spikes;
UIDs = p.Results.UIDs;
padding = p.Results.padding;
saveMat = p.Results.saveMat;
includeIntervals = p.Results.includeIntervals;

% Get session info
basename = basenameFromBasepath(basepath);
% load([basepath filesep basename '.sessionInfo.mat']);
% SR = sessionInfo.rates.wideband;

% Default events, UIDs and spikes
if isempty(spikes)
    spikes = load([basepath filesep basename '.spikes.cellinfo.mat']);
    spikes = spikes.spikes;
end
if isempty(UIDs)
    UIDs = ones(size(spikes.UID));
end
if isempty(events)
    events = load([basepath filesep basename '.ripples.events.mat']);
    events = events.ripples;
end
% Starting and ending timestamps
if isnumeric(events)
    timestamps = events;
    label = 'events';
elseif isstruct(events)
    timestamps = events.timestamps;
    label = 'events';
elseif ischar(events)
    switch lower(events)
        case 'upstates'
            UDStates = detectUD;
            events = [(UDStates.ints.DOWN(1:end-1,2)) (UDStates.ints.DOWN(2:end,2))];
            max_duration = 20;
            min_duration = 0.1;

            events(diff(events')>max_duration | diff(events')<min_duration,:) = [];
            label = 'upStates';
            
            if ~isempty(includeIntervals)
                status = InIntervals(events,includeIntervals);
                events(~status,:) = [];
            end
        case 'ripples'
            ripples = rippleMasterDetector;
            events = ripples.timestamps;
            
            if ~isempty(includeIntervals)
                status = InIntervals(events,includeIntervals);
                events(~status,:) = [];
            end
            label = 'ripples';
    end

    timestamps = events;
else
    warning('Events must be either a Nx2 vector or a bz event structure!');
end

%% Get spikes for each unit and each event

% We will save spike times of different units in different ways:
spkEventTimes = {};
% 1. Absolute and relative time of spikes by unit and by event
if size(UIDs,1) < size(UIDs,2) 
    UIDs = UIDs';
end
for unit = find(UIDs)'
    for event = 1:length(timestamps)
        % Start and end of event
        tini = timestamps(event,1) - padding;
        tend = timestamps(event,2) + padding;
        % Spikes of this unit within this event interval
        tsUnitEvent = spikes.times{unit};
        tsUnitEvent = tsUnitEvent(tsUnitEvent>=tini & tsUnitEvent<=tend);
        % Absolute time of spikes by unit and by event
        spkEventTimes.UnitEventAbs{unit,event} = tsUnitEvent';
        % Relative time of spikes by unit and by event to event start
        spkEventTimes.UnitEventRel{unit,event} = tsUnitEvent' - tini;
    end
end

% 2. Absolute and relative time of spikes by unit
for unit = find(UIDs)'
    spkEventTimes.UnitAbs{unit} = cell2mat(spkEventTimes.UnitEventAbs(unit,:));
    spkEventTimes.UnitRel{unit} = cell2mat(spkEventTimes.UnitEventRel(unit,:));
end

% 3. Absolute and relative time of spikes by event
for event = 1:length(timestamps)
    spkEventTimes.EventAbs{event} = [];
    spkEventTimes.EventRel{event} = [];
    for unit = find(UIDs)'
        spkEventTimes.EventAbs{event} = [ spkEventTimes.EventAbs{event}, ...
                                        [cell2mat(spkEventTimes.UnitEventAbs(unit,event)); ...
                                         cell2mat(spkEventTimes.UnitEventAbs(unit,event))*0+spikes.UID(unit); ...
                                         cell2mat(spkEventTimes.UnitEventAbs(unit,event))*0+unit] ];
        spkEventTimes.EventRel{event} = [ spkEventTimes.EventRel{event}, ...
                                         [cell2mat(spkEventTimes.UnitEventRel(unit,event)); ...
                                          cell2mat(spkEventTimes.UnitEventRel(unit,event))*0+spikes.UID(unit); ...
                                          cell2mat(spkEventTimes.UnitEventAbs(unit,event))*0+unit] ];
    end
    spkEventTimes.EventAbs{event} = sortrows(spkEventTimes.EventAbs{event}')';
    spkEventTimes.EventRel{event} = sortrows(spkEventTimes.EventRel{event}')';
end

% 4. Get spike position in rank and event id
for ii = 1:size(spkEventTimes.UnitEventAbs,1)
    spikeNumber = [];
    eventID = [];
    for jj = 1:size(spkEventTimes.UnitEventAbs,2)
        spikeNumber = [spikeNumber 1:length(spkEventTimes.UnitEventAbs{ii,jj})];
        eventID = [eventID ones(1,length(spkEventTimes.UnitEventAbs{ii,jj}))*jj];
    end
    spkEventTimes.spikesNumberInEvent{ii} =  spikeNumber;
    spkEventTimes.eventID{ii} =  eventID;
end

% 5. Compute stats
for ii = 1:length(spkEventTimes.UnitAbs)
    % first spike
    target_times = spkEventTimes.UnitAbs{ii}(spkEventTimes.spikesNumberInEvent{ii}==1);
    spkEventTimes.mean_firstSpike_AbsTime(ii,1) = mean(target_times);
    spkEventTimes.median_firstSpike_AbsTime(ii,1) = median(target_times);
    spkEventTimes.std_firstSpike_AbsTime(ii,1) = std(target_times);
    spkEventTimes.ci95_firstSpike_AbsTime(ii,1) = 1.96 * std(target_times)/sqrt(length(target_times));
    
    target_times = spkEventTimes.UnitRel{ii}(spkEventTimes.spikesNumberInEvent{ii}==1);
    spkEventTimes.mean_firstSpike_RelTime(ii,1) = mean(target_times);
    spkEventTimes.median_firstSpike_RelTime(ii,1) = median(target_times);
    spkEventTimes.std_firstSpike_RelTime(ii,1) = std(target_times);
    spkEventTimes.ci95_firstSpike_RelTime(ii,1) = 1.96 * std(target_times)/sqrt(length(target_times));
    
    % all spikes
    target_times = spkEventTimes.UnitAbs{ii};
    spkEventTimes.mean_allSpikes_AbsTime(ii,1) = mean(target_times);
    spkEventTimes.median_allSpikes_AbsTime(ii,1) = median(target_times);
    spkEventTimes.std_allSpikes_AbsTime(ii,1) = std(target_times);
    spkEventTimes.ci95_allSpikes_AbsTime(ii,1) = 1.96 * std(target_times)/sqrt(length(target_times));
    
    target_times = spkEventTimes.UnitRel{ii};
    spkEventTimes.mean_allSpikes_RelTime(ii,1) = mean(target_times);
    spkEventTimes.median_allSpikes_RelTime(ii,1) = median(target_times);
    spkEventTimes.std_allSpikes_RelTime(ii,1) = std(target_times);
    spkEventTimes.ci95_allSpikes_RelTime(ii,1) = 1.96 * std(target_times)/sqrt(length(target_times));
    
    % mean of event
    target_times_abs = [];
    target_times_rel = [];
    for jj = 1:size(spkEventTimes.UnitEventAbs,2)
        target_times_abs(jj) = mean(spkEventTimes.UnitEventAbs{ii,jj});
        target_times_rel(jj) = mean(spkEventTimes.UnitEventRel{ii,jj});
    end
    spkEventTimes.mean_allEventsMean_RelTime(ii,1) = nanmean(target_times_rel);
    spkEventTimes.median_allEventsMean_RelTime(ii,1) = nanmedian(target_times_rel);
    spkEventTimes.std_allEventsMean_RelTime(ii,1) = nanstd(target_times_rel);
    spkEventTimes.ci95_allEventsMean_RelTime(ii,1) = 1.96 * nanstd(target_times_rel)/sqrt(length(find(~isnan(target_times_rel))));

    spkEventTimes.mean_allEventsMean_AbsTime(ii,1) = nanmean(target_times_abs);
    spkEventTimes.median_allEventsMean_AbsTime(ii,1) = nanmedian(target_times_abs);
    spkEventTimes.std_allEventsMean_AbsTime(ii,1) = nanstd(target_times_abs);
    spkEventTimes.ci95_allEventsMean_AbsTime(ii,1) = 1.96 * nanstd(target_times_abs)/sqrt(length(find(~isnan(target_times_abs))));

    % median of event
    target_times_abs = [];
    target_times_rel = [];
    for jj = 1:size(spkEventTimes.UnitEventAbs,2)
        target_times_abs(jj) = median(spkEventTimes.UnitEventAbs{ii,jj});
        target_times_rel(jj) = median(spkEventTimes.UnitEventRel{ii,jj});
    end
    spkEventTimes.mean_allEventsMedian_RelTime(ii,1) = nanmean(target_times_rel);
    spkEventTimes.median_allEventsMedian_RelTime(ii,1) = nanmedian(target_times_rel);
    spkEventTimes.std_allEventsMedian_RelTime(ii,1) = nanstd(target_times_rel);
    spkEventTimes.ci95_allEventsMedian_RelTime(ii,1) = 1.96 * nanstd(target_times_rel)/sqrt(length(find(~isnan(target_times_rel))));

    spkEventTimes.mean_allEventsMedian_AbsTime(ii,1) = nanmean(target_times_abs);
    spkEventTimes.median_allEventsMedian_AbsTime(ii,1) = nanmedian(target_times_abs);
    spkEventTimes.std_allEventsMedian_AbsTime(ii,1) = nanstd(target_times_abs);
    spkEventTimes.ci95_allEventsMedian_AbsTime(ii,1) = 1.96 * nanstd(target_times_abs)/sqrt(length(find(~isnan(target_times_abs))));
end


% Save
if saveMat
   save([basepath filesep basename '.spikesRank_' label '.cellinfo.mat'],'spkEventTimes'); 
end


end