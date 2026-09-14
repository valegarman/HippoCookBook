function eventCSD = compute_eventCSD(lfp, events, varargin)
% COMPUTE_EVENTCSD  Event-triggered LFP average and CSD, computed per shank.
%
% INPUTS:
%   lfp    : struct with .data/.timestamps/.samplingRate, OR numeric
%            matrix [samples x channels], OR cell array of numeric
%            matrices (multiple trials).
%   events : event times (seconds)
%
% Name-Value pairs:
%   'skipChannels' : raw channel IDs to exclude and interpolate (default [])
%   'samplingRate'  : (default 1250)
%   'twin'          : [pre post] window around event, in seconds (default [0.1 0.1])
%   'spat_smooth'   : spatial smoothing window, in channels (default 11)
%   'temp_smooth'   : temporal smoothing window, in samples (default 11)
%   'doDetrend'     : detrend each channel before CSD (default false)
%   'plotCSD'       : plot CSD+LFP overlay per shank (default true)
%   'saveMat'       : save output to disk (default true)
%   'session'       : preloaded session struct (default [], loads from disk if needed)

% 2026 MV

p = inputParser;
addParameter(p,'skipChannels',[],@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_smooth',11,@isnumeric);
addParameter(p,'temp_smooth',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'session',[],@loadSession);
parse(p,varargin{:});
skipChannels = p.Results.skipChannels;
samplingRate = p.Results.samplingRate;
spat_smooth  = p.Results.spat_smooth;
temp_smooth  = p.Results.temp_smooth;
doDetrend    = p.Results.doDetrend;
doPlot       = p.Results.doPlot;
saveMat      = p.Results.saveMat;
session      = p.Results.session;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end
twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Compute event-triggered LFP average
% NOTE: always use the FULL channel set so that elecGroups channel IDs
% keep matching the columns of lfp_avg. Bad channels are handled
% afterwards, per shank, via interpolation/extrapolation.
events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
nChan = size(data,2);
lfp_temp = nan(twin(1)+twin(2)+1,nChan,length(events));
for e = 1:length(events)
    lfp_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),:);
end
lfp_avg = nanmean(lfp_temp,3)*-1;

%% Shank grouping
if isempty(session)
    if exist([basenameFromBasepath(pwd) '.session.mat']) == 2
        session = loadSession;
        elecGroups = session.extracellular.electrodeGroups.channels;
    else
        elecGroups = {1:size(lfp_avg,2)};
    end
else
    elecGroups = session.extracellular.electrodeGroups.channels;
end

% Merge user-provided skipChannels with any 'Bad' channels tagged in the
% session, if present.
if isfield(session,'channelTags') && isfield(session.channelTags,'Bad') ...
        && isfield(session.channelTags.Bad,'channels') && ~isempty(session.channelTags.Bad.channels)
    badFromSession = session.channelTags.Bad.channels;
    nBadNew = numel(setdiff(badFromSession, skipChannels));
    if nBadNew > 0
        fprintf('Adding %d bad channel(s) from session.channelTags.Bad to skipChannels.\n', nBadNew);
    end
    skipChannels = union(skipChannels, badFromSession);
end

nShanks = numel(elecGroups);
eventCSD = struct();
taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;

if doPlot
    figure('Position',[100 100 350*nShanks 500]);
end

for s = 1:nShanks
    chIDs = elecGroups{s};
    shankLFP = lfp_avg(:, chIDs);          % [time x nChanShank], top->bottom order

    isBad = ismember(chIDs, skipChannels);
    if any(isBad)
        if all(isBad)
            warning('Shank %d: all channels are in skipChannels, cannot interpolate.', s);
        else
            goodPos = find(~isBad);
            badPos  = find(isBad);
            shankLFP(:,badPos) = interp1(goodPos, shankLFP(:,goodPos)', badPos, 'linear', 'extrap')';
        end
    end

    % --- detrend ---
    if doDetrend
        shankLFP = detrend(shankLFP')';
    end

    % --- temporal smoothing ---
    if temp_smooth > 0
        for ch = 1:size(shankLFP,2)
            shankLFP(:,ch) = smooth(shankLFP(:,ch), temp_smooth, 'sgolay');
        end
    end

    % --- spatial smoothing (within this shank only) ---
    if spat_smooth > 0 && size(shankLFP,2) > 1
        for t = 1:size(shankLFP,1)
            shankLFP(t,:) = smooth(shankLFP(t,:), spat_smooth, 'lowess');
        end
    end

    % --- CSD (2nd spatial derivative, within this shank only) ---
    if size(shankLFP,2) >= 3
        shankCSD = diff(shankLFP,2,2);     % [time x (nChanShank-2)]
    else
        shankCSD = [];
        warning('Shank %d: fewer than 3 channels, cannot compute CSD.', s);
    end

    eventCSD.shankIdx(s)    = s;
    eventCSD.chanIDs{s}     = chIDs;
    eventCSD.badChannels{s} = chIDs(isBad);
    eventCSD.lfp_data{s}    = shankLFP;
    eventCSD.csd{s}         = shankCSD;
    eventCSD.nChannels(s)   = numel(chIDs);

    %% Plot: CSD as colormap, LFP traces overlaid on top
    if doPlot && ~isempty(shankCSD)
        subplot(1,nShanks,s);
        cmax = max(abs(shankCSD(:)));
        % CSD channels are offset by 1 relative to shankLFP (diff,2 drops
        % the first and last channel), so align y-axis in "channel index"
        % units consistent with shankLFP's spacing.
        contourf(taxis, 2:size(shankLFP,2)-1, shankCSD', 40, 'LineColor','none');
        hold on;
        colormap(jet); caxis([-cmax cmax]);
        set(gca,'YDir','reverse');

        % overlay LFP traces, scaled/offset to align with channel axis
        lfp_scale = 0.5 / max(abs(shankLFP(:)));  % scale traces to ~half a channel-spacing
        for ch = 1:size(shankLFP,2)
            plot(taxis, ch + lfp_scale*shankLFP(:,ch), 'k', 'LineWidth', 1);
        end
        plot([0 0], [1 size(shankLFP,2)], '--w');

        xlabel('time (ms)'); ylabel('channel (depth order)');
        title(sprintf('Shank %d', s));
        xlim([taxis(1) taxis(end)]);
        ylim([0.5 size(shankLFP,2)+0.5]);
    end
end

eventCSD.timestamps   = taxis;
eventCSD.samplingRate = samplingRate;
eventCSD.params.spat_smooth  = spat_smooth;
eventCSD.params.temp_smooth  = temp_smooth;
eventCSD.params.detrend      = doDetrend;
eventCSD.params.skipChannels = skipChannels;

if saveMat
    save([basenameFromBasepath(pwd) '.eventCSD.channelinfo.mat'],'eventCSD');
end

end