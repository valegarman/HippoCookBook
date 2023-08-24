function [psth] = spikesPsth(timestamps,varargin)
% Computes Psth for timestamps entered as inputs.
% USAGE
%   [psth] = spikesPsth(timestamps,<options>)
%
% INPUTS
%   timestamps - mx2 matrix indicating timetamps (in seconds) over which to
%                   compute psth
%   
% <OPTIONALS>
%   basepath - default pwd
%   spikes - buzcode spikes structure
%   numRep - For bootstraping, default 500. If 0, no bootstraping
%   binSize - In seconds, default 0.001 (1 ms)
%   winSize - In seconds, default 0.5 (500 ms)
%   rasterPlot - Default true
%   ratePlot - Default true
%   saveMat - default true
%   eventType - default, date, other options: slowOscillations, ripples...
%   event_ints - interval around events timestamps to compute cell reponses
%   baseline_ints - interval before events timestamps to compute baseline
%   min_pulsesNumber - minimum number of pulses to create pulses entry, default 100
%   win_Z - Interval arround events to compting Zscore mean and SD. By
%       default [-winSize/2 -events_ints(1)];
%
% OUTPUTS
%   psth - struct
%
% To do: implement post-baseline
% Developed by Pablo Abad and Manuel Valero 2022.
% Modified by Pablo Abad to restrict analysis to specific intervals
% (restrictIntervals)
%% Defaults and Params
p = inputParser;

addRequired(p,'timestamps',@isnumeric);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'numRep',500,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',1,@isnumeric);
addParameter(p,'rasterPlot',true,@islogical);
addParameter(p,'ratePlot',true,@islogical);
addParameter(p,'winSizePlot',[-.1 .5],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'eventType',date,@ischar);
addParameter(p,'event_ints',[-0.02 0.02],@isnumeric);
addParameter(p,'baseline_ints',[-0.5 -0.46],@isnumeric);
addParameter(p,'min_pulsesNumber',100,@isnumeric);
addParameter(p,'win_Z',[],@isnumeric);
addParameter(p,'win_raster',[-.5 .5],@isnumeric);
addParameter(p,'binSize_raster',[0.001],@isnumeric);
addParameter(p,'restrictIntervals',[],@isnumeric);

parse(p, timestamps,varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
rasterPlot = p.Results.rasterPlot;
ratePlot = p.Results.rasterPlot;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
savePlot = p.Results.savePlot;
force = p.Results.force;
eventType = p.Results.eventType;
event_ints = p.Results.event_ints;
baseline_ints = p.Results.baseline_ints;
min_pulsesNumber = p.Results.min_pulsesNumber;
win_Z = p.Results.win_Z;
win_raster = p.Results.win_raster;
binSize_raster = p.Results.binSize_raster;
restrictIntervals = p.Results.restrictIntervals;

%% Session Template
% Deal with inputs
prevPath = pwd;
cd(basepath);

session = loadSession;
if (exist([session.general.name '.' eventType '_psth.cellinfo.mat'],'file') || ...
        exist([session.general.name '.' eventType '_psth.cellinfo.mat'],'file')) ...
        && ~force
    disp(['Psth already computed for ', session.general.name, ' ', eventType,'. Loading file.']);
    load([session.general.name '.' eventType '_psth.cellinfo.mat']);
    return
end

if isempty(win_Z)
    win_Z = [-winSize/2 -event_ints(1)];
end

% default detection parameters
if strcmpi(eventType,'slowOscillations')
    if isempty(timestamps)
        UDStates = detectUpsDowns;
        timestamps = UDStates.timestamps.DOWN;
    end
    warning('Using default parameters for slow oscillations!');
    binSize = 0.01;
    winSize = 1;
    winSizePlot = [-0.5 0.5];
    event_ints = [-0.05 0.05];
    baseline_ints = [-0.5 0.4];
    win_Z = [-0.5 -0.1];
elseif strcmpi(eventType,'ripples')
    if isempty(timestamps)
        ripples = rippleMasterDetector;
        timestamps = ripples.peaks;
    end
    warning('Using default parameters for ripples!');
    binSize = 0.005;
    winSize = 1;
    winSizePlot = [-0.5 0.5];
    event_ints = [-0.025 0.025];
    baseline_ints = [-0.5 -0.5 + diff(event_ints)];
    win_Z = [-0.5 -0.1];
end

%% Spikes
if isempty(spikes)
    spikes = loadSpikes();
    if ~isempty(restrictIntervals)
        for ii = 1:length(spikes.UID)
            [status] = InIntervals(spikes.times{ii},restrictIntervals);
            spikes.times{ii} = spikes.times{ii}(status);
        end
    end
end

%% Get cell response
psth = [];
timestamps_recording = timestamps(1):1/1250:timestamps(end);
nConditions = size(timestamps,2);
% We can implement different conditions (if timestamps : mxn instead mx1)
% TO DO ??
disp('Generating bootstrap template...');
nEvents = int32(size(timestamps,1));
randomEvents = [];
for i = 1:numRep
    randomEvents{i} = sort(randsample(timestamps_recording,nEvents))';
end

t = [];
disp('Computing responses...');
for ii = 1:length(spikes.UID)
    fprintf(' **Events from unit %3.i/ %3.i \n',ii, size(spikes.UID,2));
    if numRep > 0 & ~isnan(timestamps) 
        [stccg, t] = CCG([spikes.times{ii} randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
        for jj = 1:nConditions
%             t_duringPulse = t > 0 & t < conditions(jj,1);
            t_duringPulse = t > event_ints(1) & t < event_ints(2);
            randomRatesDuringPulse = nanmean(stccg(t_duringPulse,2:size(randomEvents,2)+1,1),1);
            psth.bootsTrapRate(ii,jj) = mean(randomRatesDuringPulse);
            psth.bootsTrapRateStd(ii,jj) = std(randomRatesDuringPulse);
            pd = fitdist(randomRatesDuringPulse','normal');
            psth.bootsTrapCI(ii,jj,:) = pd.icdf([.001 0.999]);
        end
    else
        psth.bootsTrapRate(ii,1:nConditions) = NaN;
        psth.bootsTrapRateStd(ii,1:nConditions) = NaN;
        psth.bootsTrapCI(ii,1:nConditions,:) = nan(nConditions,2);
    end
    for jj = 1:nConditions
        nPulses = length(timestamps);
        if nPulses > min_pulsesNumber
            [stccg, t] = CCG({spikes.times{ii}, timestamps},[],'binSize',binSize,'duration',winSize,'norm','rate');
            psth.responsecurve(ii,jj,:) = stccg(:,2,1);
            psth.responsecurveSmooth(ii,jj,:) = smooth(stccg(:,2,1));
            t_duringPulse = t > event_ints(1) & t < event_ints(2);
            t_beforePulse = t > baseline_ints(1) & t < baseline_ints(2);
            t_Z = t<=-0.1;
            psth.responsecurveZ(ii,jj,:) = (stccg(:,2,1) - mean(stccg(t_Z,2,1)))/std(stccg(t_Z,2,1));
            psth.responsecurveZSmooth(ii,jj,:) = smooth((stccg(:,2,1) - mean(stccg(t_Z,2,1)))/std(stccg(t_Z,2,1)));
            psth.rateDuringPulse(ii,jj,1) = mean(stccg(t_duringPulse,2,1));
            psth.rateBeforePulse(ii,jj,1) = mean(stccg(t_beforePulse,2,1));
            psth.rateZDuringPulse(ii,jj,1) = mean(squeeze(psth.responsecurveZ(ii,jj,t_duringPulse)));
            [h, psth.modulationSignificanceLevel(ii,jj,1)] = kstest2(stccg(t_duringPulse,2,1),stccg(t_beforePulse,2,1));
            ci = squeeze(psth.bootsTrapCI(ii,jj,:));
            
            % Boostrap test
            if psth.rateDuringPulse(ii,jj,1) > ci(2)
                test = 1;
            elseif psth.rateDuringPulse(ii,jj,1) < ci(1)
                test = -1;
            else
                test = 0;
            end
            psth.bootsTrapTest(ii,jj,1) = test;
            
            % z-score change test
            if mean(psth.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif mean(psth.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                test = -1;
            else
                test = 0;
            end
            psth.zscoreTest(ii,jj,1) = test;
            
            % 3 ways test. If not boostrap, it would be 2 ways.
            if (psth.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2))) && psth.modulationSignificanceLevel(ii,jj,1)<0.05...
                    && mean(psth.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif (psth.rateDuringPulse(ii,jj,1) < ci(1) || isnan(ci(1))) && psth.modulationSignificanceLevel(ii,jj,1)<0.05 ...
                    && mean(psth.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                test = -1;
            else
                test = 0;
            end
            psth.threeWaysTest(ii,jj,1) = test;
            
            % raster plot
            rasterX = [];
            rasterY = [];
            for zz = 1:size(timestamps,1)
                temp_spk = spikes.times{ii}(find(spikes.times{ii} - timestamps(zz,1)  > win_raster(1) & spikes.times{ii} - timestamps(zz,1)  < win_raster(2))) - timestamps(zz,1);
                rasterX = [rasterX; temp_spk];
                if ~isempty(temp_spk)
                    rasterY = [rasterY; zz * ones(size((temp_spk)))];
                end
            end
            
            if ~isempty(rasterX)
                [rasterHist3,c] = hist3([rasterY rasterX],{1:size(timestamps,1) win_raster(1):binSize_raster:win_raster(2)});
                
                psth.raster.rasterCount{ii,jj} = rasterHist3;
                psth.raster.rasterProb{ii,jj} = rasterHist3/sum(rasterHist3(:));
                psth.raster.TrialsNumber{ii,jj} = c{1};
                psth.raster.times{ii,jj} = c{2};
                psth.raster.rasterTrials{ii,jj} = rasterY;
                psth.raster.rasterSpikesTimes{ii,jj} = rasterX;
            else
                psth.raster.rasterCount{ii,jj} = NaN;
                psth.raster.rasterProb{ii,jj} = NaN;
                psth.raster.TrialsNumber{ii,jj} = NaN;
                psth.raster.times{ii,jj} = NaN;
                psth.raster.rasterTrials{ii,jj} = NaN;
                psth.raster.rasterSpikesTimes{ii,jj} = NaN;
            end
            
        else
            %psth.responsecurve(ii,jj,:) = nan(duration/binSize + 1,1);
            %psth.responsecurveZ(ii,jj,:) = nan(duration/binSize + 1,1);
            psth.responsecurve(ii,jj,:) = NaN;
            psth.responsecurveZ(ii,jj,:) = NaN;
            psth.modulationSignificanceLevel(ii,jj,1) = NaN;
            psth.rateDuringPulse(ii,jj,1) = NaN;
            psth.rateBeforePulse(ii,jj,1) = NaN;
            psth.rateZDuringPulse(ii,jj,1) = NaN;
            psth.bootsTrapTest(ii,jj,1) = NaN;
            psth.zscoreTest(ii,jj,1) = NaN;
            psth.threeWaysTest(ii,jj,1) = NaN;

            psth.raster.rasterCount{ii,jj} = NaN;
            psth.raster.rasterProb{ii,jj} = NaN;
            psth.raster.TrialsNumber{ii,jj} = NaN;
            psth.raster.times{ii,jj} = NaN;
            psth.raster.rasterTrials{ii,jj} = NaN;
            psth.raster.rasterSpikesTimes{ii,jj} = NaN;
        end
    end
    psth.timestamps = t;
end

% Some metrics reponses
responseMetrics = [];
t = psth.timestamps;
if any(~isnan(psth.responsecurve))
    for ii = 1:length(spikes.UID)
        for jj = 1:nConditions
            t_duringPulse = t > 0 & t < 0.1; 
            responseMetrics.maxResponse(ii,jj) = max(squeeze(psth.responsecurve(ii,jj,t_duringPulse)));
            responseMetrics.minResponse(ii,jj) = min(squeeze(psth.responsecurve(ii,jj,t_duringPulse)));
            responseMetrics.maxResponseZ(ii,jj) = max(squeeze(psth.responsecurveZ(ii,jj,t_duringPulse)));
            responseMetrics.minResponseZ(ii,jj) = min(squeeze(psth.responsecurveZ(ii,jj,t_duringPulse)));

            responseCurveZ = squeeze(psth.responsecurveZSmooth(ii,jj,:));
            responseCurveZ(t<0) = 0;

            targetSD = -2;
            temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg2SD(ii,jj) = temp(1);
            temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg2SD(ii,jj) = temp(1) - 0.1;
            responseMetrics.responseDurationNeg2SD(ii,jj) = temp(1) - responseMetrics.latencyNeg2SD(ii,jj);

            targetSD = -1.5;
            temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg1_5SD(ii,jj) = temp(1);
            temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg1_5SD(ii,jj) = temp(1) - 0.1;
            responseMetrics.responseDurationNeg1_5SD(ii,jj) = temp(1) - responseMetrics.latencyNeg1_5SD(ii,jj);

            targetSD = -1;
            temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg1SD(ii,jj) = temp(1);
            temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg1SD(ii,jj) = temp(1) - 0.1;
            responseMetrics.responseDurationNeg1SD(ii,jj) = temp(1) - responseMetrics.latencyNeg1SD(ii,jj);

            targetSD = -.5;
            temp = [t(find(diff(responseCurveZ<targetSD) == 1)+1); NaN]; responseMetrics.latencyNeg_5SD(ii,jj) = temp(1);
            temp = [t(find(diff(responseCurveZ<targetSD) == -1)+1); NaN]; responseMetrics.recoveryNeg_5SD(ii,jj) = temp(1) - 0.1;
            responseMetrics.responseDurationNeg_5SD(ii,jj) = temp(1) - responseMetrics.latencyNeg_5SD(ii,jj);
        end
    end
end

psth.responseMetrics = responseMetrics;
psth.eventType = eventType;
psth.numRep = numRep;
psth.binSize = binSize;
psth.winSize = winSize;
psth.event_ints = event_ints;
psth.baseline_ints = baseline_ints;

% squeeze matrix, if number of conditions is 1
if nConditions == 1
    nameOfFields = fieldnames(psth);
    for ii = 1:length(nameOfFields)
        if ndims(psth.(nameOfFields{ii})) > 2
            psth.(nameOfFields{ii}) = squeeze(psth.(nameOfFields{ii}));
        end
    end
end

if saveMat
    try
        disp('Saving results...');
        save([basenameFromBasepath(pwd) '.' eventType '_psth.cellinfo.mat'],'psth');
    catch
        disp('Saving results...');
        save([basenameFromBasepath(pwd) '.' eventType '_psth.cellinfo.mat'],'psth','-v7.3');
    end
    disp('Saving results...');

    raster = psth.raster;
    psth = rmfield(psth,'raster');
    save([basenameFromBasepath(pwd) '.' eventType '_psth.cellinfo.mat'],'psth','-v7.3');
    save([basenameFromBasepath(pwd) '.' eventType '_raster.cellinfo.mat'],'raster','-v7.3');
end

% PLOTS
% 1. Rasters plot
if rasterPlot && any(any(~isnan(psth.responsecurve)))
    t = psth.timestamps;
    st = timestamps;
    if length(st) > 5000 % if more than 5000
        st = randsample(st, 5000);
        st = sort(st);
    end
    disp('   Plotting spikes raster and psth...');
    % [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);
    figure;
    set(gcf,'Position',[200 -500 2500 1200]);
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Events from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        rast_x = []; rast_y = [];
        for kk = 1:length(st)
            temp_rast = spikes.times{jj} - st(kk);
            temp_rast = temp_rast(temp_rast>winSizePlot(1) & temp_rast<winSizePlot(2));
            rast_x = [rast_x temp_rast'];
            rast_y = [rast_y kk*ones(size(temp_rast))'];
        end

        % spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
        resp = psth.responsecurveSmooth(jj,:);
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(rast_x, rast_y,'.','MarkerSize',1)
        hold on
        plot(t(t>winSizePlot(1) & t<winSizePlot(2)), resp(t>winSizePlot(1) & t<winSizePlot(2)) * kk/max(resp)/2,'k','LineWidth',2);
        xlim([winSizePlot(1) winSizePlot(2)]); ylim([0 kk]);
        if psth.threeWaysTest(jj) == 0
            title(num2str(jj),'FontWeight','normal','FontSize',10);
        elseif psth.threeWaysTest(jj) == -1
            title([num2str(jj) '(-)'],'FontWeight','normal','FontSize',10);
        elseif psth.threeWaysTest(jj) == 1
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
    if savePlot
        saveas(gcf,['SummaryFigures\spikesPsthRaster_', eventType ,'.png']); 
    end
end
% 2. Rate plot
if ratePlot && any(any(~isnan(psth.responsecurve)))
    t = psth.timestamps;
    figure
    for ii = 1:nConditions;
        subplot(nConditions,2,1 + ii * 2 - 2)
        imagesc([t(1) t(end)],[1 size(psth.responsecurve,1)],...
            psth.responsecurveSmooth); caxis([0 10]); colormap(jet);
        set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
        if ii == 1
            title('Rate [0 to 10 Hz]','FontWeight','normal','FontSize',10);
            ylabel('Cells');
        end
        if ii == nConditions
            xlabel('Time');
        else
            set(gca,'XTick',[]);
        end

        subplot(nConditions,2,2 + ii * 2 - 2)
        imagesc([t(1) t(end)],[1 size(psth.responsecurve,1)],...
            psth.responsecurveZSmooth); caxis([-3 3]); colormap(jet);
        set(gca,'TickDir','out'); xlim([winSizePlot(1) winSizePlot(2)]);
        if ii == 1
           title('Z Rate [-3 to 3 SD]','FontWeight','normal','FontSize',10);
           ylabel('Cells');
        end
        if ii == nConditions
            xlabel('Time');
        else
            set(gca,'XTick',[]);
        end
    end
end          
if savePlot
    saveas(gcf,['SummaryFigures\spikesPsthRate_',eventType,'.png']); 
end

cd(prevPath);
end

