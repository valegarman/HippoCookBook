
function [ csd ] = computeCSD (lfp, varargin)

% [ CSD ] = bz_CSD (lfp, varargin)
% Calculates the 1D approximation of current source density (CSD) from a
% linear array of LFPs

% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   or a [timestamps x channel] 2D matrix (samplingRate will
%                       be 1250)
%                   or empty, in which case runs getLFP locally.
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. If empty take all (default)
%       win         time interval to compute CSD. If empty take all (default)
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       plotCSD     true/false. Default false.
%       plotLFP     true/false. Default false.
%    =========================================================================

% OUTPUT:
%    CSD           a buzcode structure with fields csd.data,
%                                                   csd.timestamps
%                                                   csd.samplingRate
%                                                   csd.channels 
%                                                   csd.params

% Antonio FR, 7/18
% Manuel Valero (migrate to HippoCookBook, adding interpolating, badchannels, etc), 2022
% NCL-MV, Generate error if requesting first or last channel 2024

%% Parse inputs

p = inputParser;
addParameter(p,'channels',[],@isnumeric);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'win',[],@isnumeric);
addParameter(p,'spat_sm',0,@isnumeric);
addParameter(p,'temp_sm',0,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',false,@islogical);
addParameter(p,'plotLFP',false,@islogical);
addParameter(p,'session',[],@isstruct);
addParameter(p,'interpBadChannels',true,@islogical);
addParameter(p,'interpBadChannels_fast',true,@islogical);

parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
session = p.Results.session;
win = p.Results.win;
interpBadChannels_fast = p.Results.interpBadChannels_fast;
interpBadChannels = p.Results.interpBadChannels;


% session input
if isempty(session)
    try session = loadSession;
    catch
        warning('session file not found...');
    end
end

badChannels = [];
if ~isempty(session)
    if isfield(session.channelTags,'Bad')
        badChannels = session.channelTags.Bad.channels;
    end
end

% channel input
channelToRetrieve = [];
if ~isempty(channels)
    if length(channels) > 1
        channelToRetrieve = channels;
    else 
        channelToRetrieve = channels;
        % find shank with channel
        for ii = 1:length(session.extracellular.spikeGroups.channels)
            if ismember(channels,session.extracellular.spikeGroups.channels{ii})
                channels = session.extracellular.spikeGroups.channels{ii};
            end
        end
    end
elseif isempty(channels) && ~isempty(session) && ~isnumeric(lfp)
    % use best possible channel (longest with good channels)
        for ii = 1:length(session.extracellular.spikeGroups.channels)
            nchan(ii) = length(session.extracellular.spikeGroups.channels{ii}) - ...
                sum(ismember(session.extracellular.spikeGroups.channels{ii},badChannels));
        end
        [~,idx] = max(nchan);
        channels = session.extracellular.spikeGroups.channels{idx};
        channelToRetrieve = channels;
elseif isempty(channels)  && isnumeric(lfp)
    channels = 1:size(lfp,2);
    channelToRetrieve = channels;
else
    error('Analysis was not possible!');
end

if all(ismember(channelToRetrieve,badChannels))
    error('Channels to retrieve have been defined as bad channels!');
end

%lfp input
if ~isempty(lfp) && isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
else
    if isempty(lfp)
        lfp = getLFP(channels);
    elseif isstruct(lfp)
        [~,ia,~] = intersect(lfp.channels,channels,'stable');
        lfp.data = lfp.data(ia);
        lfp.channels = lfp.channels(ia);
    end
    data = lfp.data;
    timestamps = lfp.timestamps;
end

% win input
if ~isempty(win)
    inWin = InIntervals(timestamps,win);
    timestamps = timestamps(inWin);
    data = data(inWin,:);
    clear inWin
end

% interpBadChannels input
if (interpBadChannels || interpBadChannels_fast) && any(ismember(badChannels, channels))
    badChannels = badChannels(ismember(badChannels, channels));
    badChannels = ismember(channels, badChannels);
    data_interpolated = double(data);
    nChan = 1:length(channels);
    
    if ~interpBadChannels_fast
        for ii = 1:size(data,1)
            data_interpolated(ii,badChannels) = qinterp1(nChan(~badChannels),double(data(ii,~badChannels)),nChan(badChannels));
        end
    else 
        % fast interpolation first dicard channels on the extremes
        while badChannels(1) == 1 || badChannels(end) == 1
            warning('Discarting bad channels on the extremes');
            if badChannels(1) == 1
                badChannels(1) = [];
                data(:,1) = [];
                channels(1) = [];
            elseif badChannels(end) == 1
                badChannels(end) = [];
                data(:,end) = [];
                channels(end) = [];
            end
        end
        % Then, computes the average of the channels +1 and -1.
        bad_channels_pos = find(badChannels);
        for jj = 1:length(bad_channels_pos)
            edge_channels = [bad_channels_pos(jj)-1 bad_channels_pos(jj)+1];
            data_interpolated(:,bad_channels_pos) = int16(nanmean(data_interpolated(:,edge_channels),2));
        end
    end
end



%% Compute CSD

data = data*-1;

% detrend
if doDetrend
   data = detrend(data')';
end
    
% temporal smoothing
if temp_sm > 0
   for ch = 1:size(data,2) 
       data(:,ch) = smooth(double(data(:,ch)),temp_sm,'sgolay');
   end
end

% spatial smoothing
if spat_sm > 0
   for t = 1:size(data,1) 
       data(t,:) = smooth(double(data(t,:)),spat_sm,'lowess');
   end
end

% calculate CSD 
CSD = diff(data,2,2);
CSD = [nan(size(CSD,1),1) CSD nan(size(CSD,1),1);];

if any(channelToRetrieve ~= channels) && ~(channelToRetrieve == channels(1) || channelToRetrieve == channels(end))
    CSD = CSD(:,find(channels==channelToRetrieve));
    channels = channelToRetrieve;
else
    error('Not possible to retrieve top or bottom channel!');
end

% generate output structure
csd.data = CSD;
csd.timestamps = timestamps;
csd.samplingRate = samplingRate;
csd.channels = channels; 
csd.params.spat_sm = spat_sm;
csd.params.temp_sm = temp_sm;
csd.params.detrend = doDetrend;
csd.params.badChannels = badChannels; 
csd.params.interpBadChannels = interpBadChannels;
csd.params.interpBadChannels_fast = interpBadChannels_fast;

%% Plot

if plotLFP
    
    cmax = max(max(CSD)); 

    figure;
    subplot(1,2,1);
    contourf(timestamps(win(1):win(2)),1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
   
    subplot(1,2,2);
    for ch=1:size(lfp_frag,2)
        offset = 500*(ch-1);
        sh_tmp = 10e5*(lfp_frag(:,ch)) + offset;
        plot(timestamps(win(1):win(2)),sh_tmp,'k','LineWidth',1.5); hold on;
        clear sh_tmp
    end
    set(gca,'YDir','reverse','YTickLabel',[]);ylim([-500 offset+500]);xlim([timestamps(win(1)) timestamps(win(2))]);
    xlabel('time (s)');ylabel('channel');title('LFP');   
    
elseif plotCSD  
    
     cmax = max(max(CSD)); 
   
     figure;
     contourf(timestamps(win(1):win(2)),1:size(CSD,2),CSD',40,'LineColor','none');hold on;
     colormap jet; caxis([-cmax cmax]);
     set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title(CSD); 
   
end

end


    




