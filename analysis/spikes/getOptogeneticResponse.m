
function [optogeneticResponses] = getOptogeneticResponse(varargin)
% [optogeneticResponse] = getOptogeneticResponse(varargin)
%
% Computes Psth and a several statistical measures of the cell responses.
%
% <OPTIONALS>
% analogChannelsList          List of analog channels with light pulses. By default, []
% digitalChannelsList         List of digital channels with light pulses. By default,
%                       [].
% spikes            buzcode spikes structure, if not provided tries loadSpikes.
% basepath          By default pwd.
% numRep            For boostraping, default, 500. If 0, no boostraping.
% binSize           In seconds, default, 0.001.
% winSize           In seconds, default, 0.5.
% rasterPlot        Default true.
% ratePlot          Default true.
% winSizePlot       Default [-0.1 .5];
% force             Default, false.
% saveEventsFile    Default, true.
%                   
%
% OUTPUTS
% optogeneticResponse
%
% Manu-BuzsakiLab 2021
% To do: include different test for different durations stimuli
% to do: include statistical test in figures (as with uLEDs)

% Parse options
p = inputParser;
addParameter(p,'analogChannelsList',NaN);
addParameter(p,'digitalChannelsList',NaN);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',50,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',1,@isnumeric);
addParameter(p,'rasterPlot',true,@islogical);
addParameter(p,'ratePlot',true,@islogical);
addParameter(p,'winSizePlot',[-.1 .5],@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'minNumberOfPulses',200,@isnumeric);
addParameter(p,'minDuration',0.004,@isnumeric); % 4 ms
addParameter(p,'saveEventsFile',true,@islogical);
addParameter(p,'duration_round_decimal',3,@isscalar);

parse(p, varargin{:});
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;
basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
rasterPlot = p.Results.rasterPlot;
ratePlot = p.Results.rasterPlot;
saveMat = p.Results.saveMat;
winSizePlot = p.Results.winSizePlot;
force = p.Results.force;
minNumberOfPulses = p.Results.minNumberOfPulses;
minDuration = p.Results.minDuration;
saveEventsFile = p.Results.saveEventsFile;
duration_round_decimal = p.Results.duration_round_decimal;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.optogeneticResponse.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Optogenetic responses already computed! Loading file...');
    load(targetFile.name);
    return
end

%% Digital and Analog Pulses
if isnan(analogChannelsList)
    try 
        session = loadSession(basepath);
        analogChannelsList = session.analysisTags.analog_optogenetic_channels;
    catch
        warning('There is a problem with the analog channels...');
    end
end

if isnan(digitalChannelsList)
    try
        session = loadSession(basepath);
        digitalChannelsList = session.analysisTags.digital_optogenetic_channels;
    catch
        warning('There is a problem with the digital channels...');
    end
end

%%
if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

pulsesAnalog.timestamps = []; pulsesAnalog.analogChannelsListannel = [];
if strcmpi(analogChannelsList,'all')
    pulsesAnalog = getAnalogPulses;
else
    pulsesAnalog = getAnalogPulses('analogChannelsList',analogChannelsList);
end
if isempty(pulsesAnalog)
    pulsesAnalog.timestamps = []; 
    pulsesAnalog.analogChannelsList = [];
    lastAnalogChannels = 0;
else
    lastAnalogChannels = max(pulsesAnalog.analogChannelsList);
end

pulsesDigital.timestamps = []; pulsesDigital.digitalChannelsList = [];
if ~isempty(digitalChannelsList)
    digitalIn = getDigitalIn;
    for ii = 1:length(digitalChannelsList)
        pulsesDigital.timestamps = [pulsesDigital.timestamps; digitalIn.ints{digitalChannelsList(ii)}];
        pulsesDigital.digitalChannelsList = [pulsesDigital.digitalChannelsList; ones(size(digitalIn.ints{digitalChannelsList(ii)},1),1) * digitalChannelsList(ii)];
    end
end

pulses.timestamps = [pulsesAnalog.timestamps; pulsesDigital.timestamps];  % combine pulses
pulses.channel = [pulsesAnalog.analogChannelsList; pulsesDigital.digitalChannelsList + lastAnalogChannels];  % combine pulses
pulses.analogChannelsListannel = [pulsesAnalog.analogChannelsList; nan(size(pulsesDigital.digitalChannelsList))];  % 
pulses.digitalChannelsListannel = [nan(size(pulsesAnalog.analogChannelsList)); pulsesDigital.digitalChannelsList];  % 
pulses.duration = round(pulses.timestamps(:,2) - pulses.timestamps(:,1),3);  % 
pulses.isAnalog = [ones(size(pulsesAnalog.analogChannelsList)); zeros(size(pulsesDigital.digitalChannelsList))];
pulses.isDigital = [zeros(size(pulsesAnalog.analogChannelsList)); ones(size(pulsesDigital.digitalChannelsList))];

% get cell response
optogeneticResponses = [];
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
conditions(notEnoughtPulses,:) = [];
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
             indexToMerge{i}(j) = find(conditions(:,3) == max(conditions(index{i}(j),3), conditions(index{i}(j)+1,3)));
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
        

%%
spikes = loadSpikes;
% generate random events for boostraping
disp('Generating boostrap template...');
nPulses = int32(size(pulses.timestamps,1));
randomEvents = [];
for kk = 1:numRep
    randomEvents{kk} = sort(randsample(timestamps_recording, nPulses))';
end
disp('Computing responses...');
for ii = 1:length(spikes.UID)
    fprintf(' **Pulses from unit %3.i/ %3.i \n',ii, size(spikes.UID,2)); %\n
    if numRep > 0
        [stccg, t] = CCG([spikes.times{ii} randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
        for jj = 1:nConditions
            t_duringPulse = t > 0 & t < conditions(jj,1); 
            randomRatesDuringPulse = nanmean(stccg(t_duringPulse,2:size(randomEvents,2)+1,1),1);
            optogeneticResponses.bootsTrapRate(ii,jj) = mean(randomRatesDuringPulse);
            optogeneticResponses.bootsTrapRateStd(ii,jj) = std(randomRatesDuringPulse);
            pd = fitdist(randomRatesDuringPulse','normal');
            optogeneticResponses.bootsTrapCI(ii,jj,:) = pd.icdf([.001 0.999]);
        end
    else
        optogeneticResponses.bootsTrapRate(ii,1:nConditions) = NaN;
        optogeneticResponses.bootsTrapRateStd(ii,1:nConditions) = NaN;
        optogeneticResponses.bootsTrapCI(ii,1:nConditions,:) = nan(nConditions,2);
    end
    for jj = 1:nConditions
        pul = pulses.timestamps(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1),:);
        isAnalog = median(pulses.isAnalog(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1)));
        channelPulse = median(pulses.channel(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1)));
        
        nPulses = length(pul);
        pulseDuration = conditions(jj,1);
        if nPulses > 100
            [stccg, t] = CCG({spikes.times{ii}, pul(:,1)},[],'binSize',binSize,'duration',winSize,'norm','rate');
            optogeneticResponses.responsecurve(ii,jj,:) = stccg(:,2,1);
            optogeneticResponses.responsecurveSmooth(ii,jj,:) = smooth(stccg(:,2,1));
            t_duringPulse = t > 0 & t < pulseDuration; 
            t_beforePulse = t > -pulseDuration & t < 0; 
            optogeneticResponses.responsecurveZ(ii,jj,:) = (stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1));
            optogeneticResponses.responsecurveZSmooth(ii,jj,:) = smooth((stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1)));
            optogeneticResponses.rateDuringPulse(ii,jj,1) = mean(stccg(t_duringPulse,2,1));
            optogeneticResponses.rateBeforePulse(ii,jj,1) = mean(stccg(t_beforePulse,2,1));
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = mean(squeeze(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)));
            optogeneticResponses.durationPerPulse(ii,jj,1) = t(find(t_duringPulse,1,'last')+1) - t(find(t_duringPulse,1,'first')-1);
            optogeneticResponses.isAnalog(ii,jj,1) = isAnalog;
            optogeneticResponses.isDigital(ii,jj,1) = ~isAnalog;
            optogeneticResponses.channelPulse(ii,jj,1) = channelPulse;
            
            [h, optogeneticResponses.modulationSignificanceLevel(ii,jj,1)] = kstest2(stccg(t_duringPulse,2,1),stccg(t_beforePulse,2,1));
            ci = squeeze(optogeneticResponses.bootsTrapCI(ii,jj,:));
            
            % Boostrap test
            if optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2)
                test = 1;
            elseif optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(1)
                test = -1;
            else
                test = 0;
            end
            optogeneticResponses.bootsTrapTest(ii,jj,1) = test;
            
            % z-score change test
            if mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                test = -1;
            else
                test = 0;
            end
            optogeneticResponses.zscoreTest(ii,jj,1) = test;
            
            % 3 ways test. If not boostrap, it would be 2 ways.
            if (optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2))) && optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01...
                    && mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif (optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(1) || isnan(ci(1))) && optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01 ...
                    && mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                test = -1;
            else
                test = 0;
            end
            optogeneticResponses.threeWaysTest(ii,jj,1) = test;
        else
            optogeneticResponses.responsecurve(ii,jj,:) = nan; %nan(pulseDuration/binSize + 1,1);
            optogeneticResponses.responsecurveZ(ii,jj,:) = nan; %(pulseDuration/binSize + 1,1);
            optogeneticResponses.modulationSignificanceLevel(ii,jj,1) = NaN;
            optogeneticResponses.rateDuringPulse(ii,jj,1) = NaN;
            optogeneticResponses.rateBeforePulse(ii,jj,1) = NaN;
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = NaN;
            optogeneticResponses.bootsTrapTest(ii,jj,1) = NaN;
            optogeneticResponses.zscoreTest(ii,jj,1) = NaN;
            optogeneticResponses.threeWaysTest(ii,jj,1) = NaN;
            optogeneticResponses.durationPerPulse(ii,jj,1) = NaN;
            optogeneticResponses.isAnalog(ii,jj,1) = isAnalog;
            optogeneticResponses.isDigital(ii,jj,1) = ~isAnalog;
            optogeneticResponses.channelPulse(ii,jj,1) = channelPulse;
        end
    end
    optogeneticResponses.timestamps = t;
end

% find intervals
lag = 20; % max interval between pulses, in seconds
stimulationEpochs(1,1) = pulses.timestamps(1,1);
intPeaks =find(diff(pulses.timestamps(:,1))>lag);
for ii = 1:length(intPeaks)
    stimulationEpochs(ii,2) = pulses.timestamps(intPeaks(ii),2);
    stimulationEpochs(ii+1,1) = pulses.timestamps(intPeaks(ii)+1,1);
end
stimulationEpochs(end,2) = pulses.timestamps(end,2);

optogeneticResponses.channels = conditions(:,2);
optogeneticResponses.pulseDuration = conditions(:,1);
optogeneticResponses.pulses = pulses;
optogeneticResponses.numRep = numRep;
optogeneticResponses.binSize = binSize;
optogeneticResponses.winSize = winSize;
optogeneticResponses.winSizePlot = winSizePlot;
optogeneticResponses.conditions = conditions;
optogeneticResponses.conditionsLabels = {'durations','channels','numberOfPulses'};
optogeneticResponses.nConditions = nConditions;
optogeneticResponses.stimulationEpochs = stimulationEpochs;

% Some metrics reponses
responseMetrics = [];
t = optogeneticResponses.timestamps;
for ii = 1:length(spikes.UID)
    for jj = 1:nConditions
        t_duringPulse = t > 0 & t < optogeneticResponses.pulseDuration(jj); 
        responseMetrics.maxResponse(ii,jj) = max(squeeze(optogeneticResponses.responsecurve(ii,jj,t_duringPulse)));
        responseMetrics.minResponse(ii,jj) = min(squeeze(optogeneticResponses.responsecurve(ii,jj,t_duringPulse)));
        responseMetrics.maxResponseZ(ii,jj) = max(squeeze(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)));
        responseMetrics.minResponseZ(ii,jj) = min(squeeze(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)));
        
        responseCurveZ = squeeze(optogeneticResponses.responsecurveZSmooth(ii,jj,:));
        responseCurveZ(t<0) = 0;
        
        targetSD = -2;
        temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg2SD(ii,jj) = temp(1);
        temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg2SD(ii,jj) = temp(1) - optogeneticResponses.pulseDuration(jj);
        responseMetrics.responseDurationNeg2SD(ii,jj) = temp(1) - responseMetrics.latencyNeg2SD(ii,jj);
        
        targetSD = -1.5;
        temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg1_5SD(ii,jj) = temp(1);
        temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg1_5SD(ii,jj) = temp(1) - optogeneticResponses.pulseDuration(jj);
        responseMetrics.responseDurationNeg1_5SD(ii,jj) = temp(1) - responseMetrics.latencyNeg1_5SD(ii,jj);
        
        targetSD = -1;
        temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg1SD(ii,jj) = temp(1);
        temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg1SD(ii,jj) = temp(1) - optogeneticResponses.pulseDuration(jj);
        responseMetrics.responseDurationNeg1SD(ii,jj) = temp(1) - responseMetrics.latencyNeg1SD(ii,jj);
        
        targetSD = -.5;
        temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg_5SD(ii,jj) = temp(1);
        temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg_5SD(ii,jj) = temp(1) - optogeneticResponses.pulseDuration(jj);
        responseMetrics.responseDurationNeg_5SD(ii,jj) = temp(1) - responseMetrics.latencyNeg_5SD(ii,jj);
    end
end

optogeneticResponses.responseMetrics = responseMetrics;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.optogeneticResponse.cellinfo.mat'],'optogeneticResponses');
end

if saveEventsFile
    disp('Saving timestamps...');
    filename = split(pwd,filesep); filename = filename{end};
    optoPulses = optogeneticResponses.pulses;
    optoPulses.stimulationEpochs = stimulationEpochs;
    save([filename '.optogeneticPulses.events.mat'],'optoPulses');
end

% PLOTS
% 1. Rasters plot
if rasterPlot
    t = optogeneticResponses.timestamps;
    for ii = 1:optogeneticResponses.nConditions
        st = pulses.timestamps(pulses.channel == conditions(ii,2) & pulses.duration == conditions(ii,1),1);
        
        if length(st) > 5000 % if more than 5000
            st = randsample(st, 5000);
            st = sort(st);
        end
        disp('   Plotting spikes raster and psth...');
        % [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);
        figure;
        set(gcf,'Position',[200 -500 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            rast_x = []; rast_y = [];
            for kk = 1:length(st)
                temp_rast = spikes.times{jj} - st(kk);
                temp_rast = temp_rast(temp_rast>winSizePlot(1) & temp_rast<winSizePlot(2));
                rast_x = [rast_x temp_rast'];
                rast_y = [rast_y kk*ones(size(temp_rast))'];
            end

            % spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
            resp = squeeze(optogeneticResponses.responsecurveSmooth(jj,ii,:));
            dur = optogeneticResponses.durationPerPulse(jj,ii,1);
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            plot(rast_x, rast_y,'.','MarkerSize',1,'color',[.6 .6 .6]);
            hold on
            plot(t(t>winSizePlot(1) & t<winSizePlot(2)), resp(t>winSizePlot(1) & t<winSizePlot(2)) * kk/max(resp)/2,'k','LineWidth',2);
            xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 kk*1.1]);
            plot([0 dur],[kk*1.05 kk*1.05],'color',[0 0.6 0.6],'LineWidth',2);
            
            if optogeneticResponses.threeWaysTest(jj,ii) == 0
                title(num2str(jj),'FontWeight','normal','FontSize',10);
            elseif optogeneticResponses.threeWaysTest(jj,ii) == -1
                title([num2str(jj) '(-)'],'FontWeight','normal','FontSize',10);
            elseif optogeneticResponses.threeWaysTest(jj,ii) == 1
                title([num2str(jj) '(+)'],'FontWeight','normal','FontSize',10);
            end

            if jj == 1
                ylabel('Trial');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\OptogenticRespRaster_ch',num2str(conditions(ii,2)),'_dur',num2str(conditions(ii,1)),'ch.png']); 
        close(gcf);
    end
end
% 2. Rate plot
if ratePlot
    t = optogeneticResponses.timestamps;
    figure
    set(gcf,'Position',[200 -100 1000 700]);
    for ii = 1:nConditions
        subplot(nConditions,2,1 + ii * 2 - 2)
        imagesc([t(1) t(end)],[1 size(optogeneticResponses.responsecurve,1)],...
            squeeze(optogeneticResponses.responsecurveSmooth(:,ii,:))); caxis([0 10]); colormap(jet);
        set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
        if ii == 1
            title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
        end
        if ii == nConditions
            xlabel('Time');
        else
            set(gca,'XTick',[]);
        end
        ylabel(['ch: ' num2str(conditions(ii,2)) ', dur: ' num2str(conditions(ii,1)),'s']);
        
        ylim([-3 size(optogeneticResponses.responsecurveSmooth(:,ii,:),1)]);
        hold on
        plot([0 median(optogeneticResponses.durationPerPulse(:,ii,:))],[-1.5 -1.5],'color',[0 0.6 0.6],'LineWidth',2);
        
        subplot(nConditions,2,2 + ii * 2 - 2)
        imagesc([t(1) t(end)],[1 size(optogeneticResponses.responsecurve,1)],...
            squeeze(optogeneticResponses.responsecurveZSmooth(:,ii,:))); caxis([-3 3]); colormap(jet);
        set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
        if ii == 1
           title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
        end
        if ii == nConditions
            xlabel('Time');
        else
            set(gca,'XTick',[]);
        end
        ylabel(['ch: ' num2str(conditions(ii,2)) ', dur: ' num2str(conditions(ii,1)),'s']);
        
        ylim([-3 size(optogeneticResponses.responsecurveSmooth(:,ii,:),1)]);
        hold on
        plot([0 median(optogeneticResponses.durationPerPulse(:,ii,:))],[-1.5 -1.5],'color',[0 0.6 0.6],'LineWidth',2);
    end
    saveas(gcf,['SummaryFigures\optogeneticPulsesPsth.png']); 
end          


cd(prevPath);
end