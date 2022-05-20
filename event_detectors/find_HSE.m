function HSE = find_HSE(basepath,spikes,varargin)
% Find high-sychrony events among spikes of a set of units (e.g. for
% decoding).
% For each session the combined spiking of all recorded 
% CA1 cells were binned in 1ms bins and convolved with a 15 ms Gaussian 
% kernel (5). For each session a trigger rate was defined as being 3 standard deviations 
% above the mean of all 1 ms bins within NREM epochs of both PRE and POST epochs 
% combined (Grosmark 2016)
%
% INPUTS
%   'spikes'    buzcode compatible 'spikes.cellinfo' struct
%
%   (optional)
%	'algorithm'	Currently supported: 'bayes', 'PVcorr', default 'bayes'
%
% OUTPUT
%   Outvar:   	Description
%
% EXAMPLE CALLS
% [] = DefaultTemplateBuzcode(pwd)
%
% Thomas Hainmueller, 2020, Buzsakilab
% Winnie Yang, 2022
%% Input handling
p = inputParser;
addParameter(p,'nSigma',2.5,@isnumeric);
addParameter(p,'tSmooth',.015,@isnumeric); % in s
addParameter(p,'binsz',.001,@isnumeric); % in s
addParameter(p,'cell_ID',[],@isnumeric);
addParameter(p,'name',[],@ischar);
addParameter(p,'basename',[],@ischar);
addParameter(p,'save_evts',true,@islogical);
addParameter(p,'mindur',.01,@isnumeric);
addParameter(p,'maxdur',1,@isnumeric);
addParameter(p,'winSize',0.5,@isnumeric);

parse(p,varargin{:})

nSigma = p.Results.nSigma;
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
cell_ID = p.Results.cell_ID;
name = p.Results.name;
basename = p.Results.basename;
save_evts = p.Results.save_evts;
mindur = p.Results.mindur;
maxdur = p.Results.maxdur;
winSize = p.Results.winSize;

%% Set defaults

save_folder = [basepath, '\', 'rippleHSE'];
if isempty(cell_ID)
    cell_ID = 1:length(spikes.UID);
end

if isempty(name)
    name = 'HSE';
end

if isempty(basename)
    basename = basenameFromBasepath(pwd);
end

%% Get spike rate over time
% winnie changed 
allspk = cat(1,spikes.times{cell_ID});
allspk = sort(allspk);
ts = 0:binsz:allspk(end);
spkhist = hist(allspk,ts);
spkmean = mean(spkhist);

fsize = tSmooth/binsz;
gfilt = fspecial('gauss',[10*fsize 1],fsize);

spkhist = conv(spkhist,gfilt,'same');
spkhist = zscore(spkhist);

evtidx = spkhist>nSigma;
evtidx = find(diff(evtidx)==1)+1;% diff helps to find the boundary of a chunk of spkhist > nSigma
belowm = spkhist<spkmean; % Logical to run faster
[startID, stopID, evtstart, evtstop, evtdur, evtamp, evtpeak] = deal(zeros(1,length(evtidx)));
%startID = 1; % Initialize to 1 for the first comparison to work



%%
for e = 1:length(evtidx)
    startID(e) = max([1 find(belowm(1:evtidx(e)),1,'last')]);
    if startID(e)>max(stopID)
        stopID(e) = min([length(belowm) evtidx(e)+find(belowm(evtidx(e):end),1,'first')]);
        evtstart(e) = startID(e)*binsz - binsz;
        evtstop(e) = stopID(e)*binsz - binsz;
        evtdur(e) = (stopID(e) - startID(e))*binsz;
        
        % Get amplitude and peak
        [amp, peakID] = max(spkhist(startID(e):stopID(e)));
        evtamp(e) = amp;
        peakID = peakID + startID(e);
        evtpeak(e) = peakID*binsz - binsz;
    end
end
    


%%
goodHSE = find(evtdur>mindur & evtdur<maxdur); % Add mindur here, if desired
evtstart = evtstart(goodHSE);
evtstop = evtstop(goodHSE);
evtdur = evtdur(goodHSE);
evtpeak = evtpeak(goodHSE);

event_meanResponse = nan(length(evtstart),1);
event_Zresponse = nan(length(evtstart),2*winSize/binsz);
[~, spks_event] = InIntervals(spikes.spindices(:,1),[evtpeak'-winSize evtpeak'+winSize]);

for e = 1:length(evtstart)
    event_spk = spikes.spindices(find(spks_event==e),1);
    
    

    fsize = tSmooth/binsz;
    gfilt = fspecial('gauss',[10*fsize 1],fsize);
    
    %spikes of events within the window 
    ts_wd = evtpeak(e)-winSize:binsz:evtpeak(e)+winSize;
    spkhist_wd = histcounts(event_spk,ts_wd);
    spkhist_wd = conv(spkhist_wd,gfilt,'same');
    spkhist_wd = zscore(spkhist_wd);
    
    %spikes of events within the event start to event end
    ts = evtstart(e):binsz:evtstop(e);
    spkhist = histcounts(event_spk,ts);
    spkhist = conv(spkhist,gfilt,'same');
    spkhist = zscore(spkhist);   
    spkmean = mean(spkhist);
    
    % event z scored response within in a 0.5s window (for plot PETH later) 
    event_Zresponse(e,:) = spkhist_wd(1:2*winSize/binsz);
    % mean event response from start to end of the event 
    event_meanResponse(e) = spkmean;
end
    


%% Remove overlapping evts.!

%% Create buzcode event structure and save it
HSE.timestamps = cat(2,evtstart',evtstop');
HSE.peaks = evtpeak';
HSE.amplitudes = evtamp;
HSE.amplitudeUnits = 'spikes';
HSE.eventID = ones(size(evtpeak));
HSE.eventIDlabels = repmat({name},length(evtpeak),1);
HSE.eventIDbinary = false(length(evtpeak),1);
HSE.duration = evtdur;
HSE.center = evtstart + evtdur/2;
HSE.eventZresponse = event_Zresponse;
HSE.eventMeanResponse = event_meanResponse;

HSE.detectorinfo.detectorname = 'find_HSE';
HSE.detectorinfo.detectionparms = [];
HSE.detectorinfo.detectionintervals = [0 Inf];
HSE.detectorinfo.detectiondate = datetime('today');
HSE.detectorinfo.winSize = winSize;

if save_evts
    cd(save_folder);
    save([save_folder, '/',basename '.' name '.mat'],'HSE');
end

%% Create FMA .evt structure and save it
% .evt (FMA standard)
if save_evts
    n = length(evtstart);
    d1 = cat(1,evtstart,evtpeak,evtstop);%DS1triad(:,1:3)';
    events1.time = d1(:);
    for i = 1:3:3*n
        events1.description{i,1} = [name ' start'];
        events1.description{i+1,1} = [name ' peak'];
        events1.description{i+2,1} = [name ' stop'];
    end
    
    SaveEvents([save_folder, '\', basename '_' name '.HSE.evt'],events1);
end


end