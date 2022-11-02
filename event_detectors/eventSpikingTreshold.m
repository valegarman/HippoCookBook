
function [events] = eventSpikingTreshold(events,spikes,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'events'           Buzcode format events (i.e. ripples) structure.
%
% <optional>
%    'basepath'         Default 'pwd'
%    'spikes'           Buzcode format spikes structure. If not provided runs loadSpikes.      
%    'events'           Structure containing the statistical test results.
%    'spikingThreshold' .5 
%    'winSize'          .5
%    'eventSize'        .01
%    'figOpt'           Default true
%    'shanksID'         Shanks ID for loading spikes from           
% 
% OUTPUS
%    'events'           Buzcode format events (i.e. ripples) structure
%                           after event/spiking thresholing 
%
% Manu Valero - BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'spikes',[]);
addParameter(p,'spikingThreshold',.5);
addParameter(p,'winSize',.5);
addParameter(p,'eventSize',.01);
addParameter(p,'figOpt',true,@islogical);
addParameter(p,'shanksID','all');

parse(p,varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
spikingThreshold = p.Results.spikingThreshold;
winSize = p.Results.winSize;
eventSize = p.Results.eventSize;
figOpt = p.Results.figOpt;
shanksID = p.Results.shanksID;

prevPath = pwd;
cd(basepath);

% 
if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if ischar(shanksID) && strcmpi(shanksID,'all')
    inShank = ones(size(spikes.times));
else
    try
        inShank = ismember(spikes.shankID, shanksID);
        fprintf(' Discarting %3.i units from %3.i session for response  validation... \n', length(find(inShank==0)), length(inShank)); %\n
    catch
        warning('Discating spikes by shank was not possible! Using all detected spikes...');
        inShank = ones(size(spikes.times));
    end
end

[spikemat] = bz_SpktToSpkmat(spikes, 'dt',0.01,'overlap',6);
spikemat.data = spikemat.data(:,inShank);
spikemat.UID = spikemat.UID(inShank);
sSpkMat = zscore(sum(spikemat.data,2)/size(spikemat.data,2));
% sSpkMat = mean(zscore(spikemat.data,[],1),2);
clear eventPopResponse

tic
toCheck = find(events.peaks+winSize<spikemat.timestamps(end));
sizeResponse = length(int32(1:winSize*2/(mean(diff(spikemat.timestamps)))-1));
eventPopResponse = zeros(length(toCheck),sizeResponse);
parfor ii = 1: length(toCheck)
    temp = sSpkMat(spikemat.timestamps>=events.peaks(ii)-winSize ...
        & spikemat.timestamps<=events.peaks(ii)+winSize);
    if length(temp) < sizeResponse
        eventPopResponse(ii,:) = zeros(1,sizeResponse);
    else
        eventPopResponse(ii,:) = temp(int32(1:sizeResponse));
    end
end
toc

t_event = linspace(-winSize,winSize,size(eventPopResponse,2));
eventResponse = zeros(size(events.peaks));
eventResponse(toCheck) = mean(eventPopResponse(:,t_event>-eventSize & t_event<eventSize),2);
[~,idx] = sort(eventResponse(toCheck));

% 
validEvents = find(eventResponse>spikingThreshold);
events.timestamps = events.timestamps(validEvents,:);
try 
    events.peaks = events.peaks(validEvents,:);
    events.peakNormedPower = events.peakNormedPower(validEvents,:);
end
events.eventSpikingParameters.spikingThreshold = spikingThreshold;
events.eventSpikingParameters.winSize = winSize;
events.eventSpikingParameters.eventSize = eventSize;

fprintf('Keeping %4.0f of %4.0f events \n',length(validEvents),length(eventResponse));


if figOpt
    figure
    subplot(1,3,[1 2])
    hold on
    imagesc(t_event,1:size(eventPopResponse,1),eventPopResponse(idx,:),[-3 3]); colormap(jet);
    plot([t_event([1 end])], length(find(eventResponse<spikingThreshold))* ones(2,1) ,'r','LineWidth',1); axis tight;
    ylabel('Events'); xlabel('Time (s)'); set(gca,'YDir','normal','TickDir','out');

    subplot(1,3,3)
    hold on
    plot(eventResponse(idx),1:size(eventPopResponse,1)); 
    plot([spikingThreshold spikingThreshold], [1 size(eventPopResponse,1)],'r');
    xlabel('Response (SD)'); ylim([1 size(eventPopResponse,1)]); set(gca,'YDir','normal','TickDir','out');
    try 
        saveas(gcf,'SummaryFigures\eventSpikingThreshold.png');
    end
end

cd(prevPath);
end