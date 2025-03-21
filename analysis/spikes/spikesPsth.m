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
%   minNumberOfPulses - minimum number of pulses to create pulses entry, default 100
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
addParameter(p,'getRaster',true,@islogical);
addParameter(p,'ratePlot',true,@islogical);
addParameter(p,'winSizePlot',[-.1 .5],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'eventType',date,@ischar);
addParameter(p,'event_ints',[-0.02 0.02],@isnumeric);
addParameter(p,'baseline_ints',[-0.5 -0.46],@isnumeric);
addParameter(p,'minNumberOfPulses',100,@isnumeric);
% addParameter(p,'win_Z',[],@isnumeric);
addParameter(p,'restrictIntervals',[],@isnumeric);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'raster_time',[-0.250 0.250],@isnumeric);
addParameter(p,'salt_win',[0.01],@isscalar);
addParameter(p,'salt_binSize',[0.001],@isscalar);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','_psth',@ischar);
addParameter(p,'save_raster_as','_raster',@ischar);
addParameter(p,'sr',[],@isnumeric);


parse(p, timestamps,varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
getRaster = p.Results.getRaster;
ratePlot = p.Results.ratePlot;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
savePlot = p.Results.savePlot;
force = p.Results.force;
eventType = p.Results.eventType;
event_ints = p.Results.event_ints;
baseline_ints = p.Results.baseline_ints;
minNumberOfPulses = p.Results.minNumberOfPulses;
% win_Z = p.Results.win_Z;
restrictIntervals = p.Results.restrictIntervals;
bootsTrapCI = p.Results.bootsTrapCI;
salt_baseline = p.Results.salt_baseline;
raster_time = p.Results.raster_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;
save_raster_as = p.Results.save_raster_as;
sr = p.Results.sr;

%% Session Template
% Deal with inputs
prevPath = pwd;
cd(basepath);

if minNumberOfPulses < 2
    error('Number of pulses should be lager than 1');
end

session = loadSession;
if exist([session.general.name '.' eventType '_psth.cellinfo.mat'],'file') ...
        && ~force
    disp(['Psth already computed for ', session.general.name, ' ', eventType,'. Loading file.']);
    load([session.general.name '.' eventType '_psth.cellinfo.mat']);
    return
end

% if isempty(win_Z)
%     win_Z = [-winSize/2 -event_ints(1)];
% end

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
    % win_Z = [-0.5 -0.1];
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
    % win_Z = [-0.5 -0.1];
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

if ~isfield(spikes, 'numcells')
    spikes.numcells = length(spikes.times); 
end

ints = [];
session = loadSession;
if isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = '_post';
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_baseline
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [0 session.epochs{ii}.startTime];
            warning('Epoch with manipulations found! Restricting analysis to baseline interval!');
        end
    end
    if isempty(ints)
        ints = [0 Inf];
    end
else
    ints = [0 Inf];
end

restrict_ints = IntersectIntervals([ints; restrict_to]);
timestamps = Restrict(timestamps, restrict_ints);

if isempty(sr)
    session = loadSession;
    if isfield(session,'extracellular') && isfield(session.extracellular,'sr')
        sr = session.extracellular.sr;
    else
        warning('Sampling rate not provided! Using 30000 as default parameter...');
        sr = 30000;
    end
end

%% Get cell response
psth = [];
if isempty(timestamps)
    return
end

timestamps_recording = timestamps(1):1/1250:timestamps(end);
nConditions = size(timestamps,2);
% We can implement different conditions (if timestamps : mxn instead mx1)
% TO DO ??

t = [];
disp('Computing responses...');
jj = 1;

% disp('Generating bootstrap template...');
nEvents = int32(size(timestamps,1));
randomEvents = [];

for i = 1:numRep
    randomEvents{i} = sort(randsample(timestamps_recording,nEvents))';
end
pulseDuration = abs(diff(event_ints));
if ~isnan(timestamps) & numRep > 0
    [stccg, t] = CCG([spikes.times randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate','Fs',1/sr);
    fprintf('\n'); %
    t_duringPulse = t > event_ints(1) & t < event_ints(2);
    randomRatesDuringPulse = squeeze(mean(stccg(t_duringPulse, length(spikes.UID)+1:end,1:length(spikes.UID)),1));
    psth.bootsTrapRate(:,jj) = mean(randomRatesDuringPulse,1);
    psth.bootsTrapRateStd(:,jj) = std(randomRatesDuringPulse,[],1);
    psth.bootsTrapRateSEM(:,jj) = std(randomRatesDuringPulse,[],1)/sqrt(numRep);
    if ~isempty(randomRatesDuringPulse)
        for ii = 1:size(randomRatesDuringPulse,2)
            pd = fitdist(randomRatesDuringPulse(:,ii),'normal');
            psth.bootsTrapCI(ii,jj,1:2) = pd.icdf(bootsTrapCI);
        end
    else
        psth.bootsTrapCI(1:size(randomRatesDuringPulse,2),1:2) = NaN;
    end
else
    psth.bootsTrapRate(1:spikes.numcells,jj) = NaN;
    psth.bootsTrapRateStd(1:spikes.numcells,jj) = NaN;
    psth.bootsTrapRateSEM(1:spikes.numcells,jj) = NaN;
    psth.bootsTrapCI(1:spikes.numcells,jj,1:2) = NaN;
end

if ~isempty(timestamps) && all(~isnan(timestamps))
    pul = timestamps;
else
    pul = [0];
end
disp('Computing responses...');
times = spikes.times; times{length(times)+1} = pul;
[stccg, t] = CCG(times,[],'binSize',binSize,'duration',winSize,'norm','rate','Fs',1/sr); fprintf('\n'); %

psth.responsecurve(:,jj,:) = squeeze(stccg(:, end , 1:end-1))';
if length(times{end}) < minNumberOfPulses
    psth.responsecurve(:,jj,:) = psth.responsecurve(:,jj,:) * NaN;
end
t_duringPulse = t > event_ints(1) & t < event_ints(2);
t_beforePulse = t > baseline_ints(1) & t < baseline_ints(2);
t_Z = t<=event_ints(1);

numberOfPulses = size(pul,1);
for ii = 1:size(psth.responsecurve,1)
    if numberOfPulses > minNumberOfPulses
        psth.responsecurveSmooth(ii,jj,:) = smooth(psth.responsecurve(ii,jj,:));
        psth.responsecurveZ(ii,jj,:) = (psth.responsecurve(ii,jj,:)...
            - mean(psth.responsecurve(ii,jj,t_Z)))...
            /std(psth.responsecurve(ii,jj,t_Z));
        psth.responsecurveZSmooth(ii,jj,:) = smooth(psth.responsecurveZ(ii,jj,:));
        psth.rateDuringPulse(ii,jj,1) = mean(psth.responsecurve(ii,jj,t_duringPulse));
        psth.rateBeforePulse(ii,jj,1) = mean(psth.responsecurve(ii,jj,t_beforePulse));
        psth.rateZDuringPulse(ii,jj,1) = mean(psth.responsecurveZ(ii,jj,t_duringPulse));
        psth.rateZBeforePulse(ii,jj,1) = mean(psth.responsecurveZ(ii,jj,t_beforePulse));
        psth.durationPerPulse(ii,jj,1) = t(find(t_duringPulse,1,'last')+1) - t(find(t_duringPulse,1,'first')-1);
        psth.pulseDuration(ii,jj,1) = pulseDuration;
        psth.condition(ii,jj,1) = jj;

        try
            [h, psth.modulationSignificanceLevel(ii,jj,1)] = ...
                 kstest2(squeeze(psth.responsecurve(ii,jj,t_duringPulse))...
                ,squeeze(psth.responsecurve(ii,jj,t_beforePulse)));
        catch
             psth.modulationSignificanceLevel(ii,jj,1) = NaN;
        end
        
        % Boostrap test
        ci = squeeze(psth.bootsTrapCI(ii,jj,:));
        if psth.rateDuringPulse(ii,jj,1) > ci(2)
            test = 1;
        elseif psth.rateDuringPulse(ii,jj,1) < ci(1)
            test = -1;
        elseif isnan(psth.rateDuringPulse(ii,jj,1))
                test = NaN;
        else
            test = 0;
        end
        psth.bootsTrapTest(ii,jj,1) = test;
        
        % z-score change test
        if mean(psth.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
            test = 1;
        elseif mean(psth.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
             test = -1;
        elseif isnan(psth.rateZDuringPulse(ii,jj,1))
            test = NaN;
        else
            test = 0;
        end
        psth.zscoreTest(ii,jj,1) = test;

        % Generating raster
        rasterX = [];
        rasterY = [];
        if getRaster
            for zz = 1:size(pul,1)
                temp_spk = spikes.times{ii}(find(spikes.times{ii} - pul(zz,1)  > raster_time(1) & spikes.times{ii} - pul(zz,1)  < raster_time(2))) - pul(zz,1);
                rasterX = [rasterX; temp_spk];
                if ~isempty(temp_spk)
                    rasterY = [rasterY; zz * ones(size((temp_spk)))];
                end
            end
        end
        if ~isempty(rasterX)
            [rasterHist3,c] = hist3([rasterY rasterX],{1:size(pul,1) raster_time(1):salt_binSize:raster_time(2)});
            psth.raster.rasterCount{ii,jj} = rasterHist3;
            psth.raster.rasterProb{ii,jj} = rasterHist3/sum(rasterHist3(:));
            psth.raster.TrialsNumber{ii,jj} = c{1};
            psth.raster.times{ii,jj} = c{2};
            psth.raster.rasterTrials{ii,jj} = rasterY;
            psth.raster.rasterSpikesTimes{ii,jj} = rasterX;

            time  = c{2};
            % running salt
            baseidx = dsearchn(time', salt_baseline');
            tidx = dsearchn(time', [0; pulseDuration*2]);      

            st = length(baseidx(1):baseidx(2));
            nmbn = round(salt_win/salt_binSize);
            v = 1:nmbn:st;
            if any((v + nmbn - 1) > st)
                error('reduce window size or baseline duration')
            end

            [psth.salt.p_value(ii,jj,1), psth.salt.I_statistics(ii,jj,1)] = salt(rasterHist3(:,baseidx(1):baseidx(2)),rasterHist3(:,tidx(1):tidx(2)),salt_binSize, salt_win);
        else
            psth.raster.rasterCount{ii,jj} = NaN;
            psth.raster.rasterProb{ii,jj} = NaN;
            psth.raster.TrialsNumber{ii,jj} = NaN;
            psth.raster.times{ii,jj} = NaN;
            psth.raster.rasterTrials{ii,jj} = NaN;
            psth.raster.rasterSpikesTimes{ii,jj} = NaN;
            psth.salt.p_value(ii,jj,1) = NaN;
            psth.salt.I_statistics(ii,jj,1) = NaN;
        end
        
        % multiple test test. If not boostrap, it would be 2 ways.
        if (psth.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2))) && psth.modulationSignificanceLevel(ii,jj,1)<0.01...
                && mean(psth.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
            test = 1;
        elseif (psth.rateDuringPulse(ii,jj,1) < ci(1) || isnan(ci(1))) && psth.modulationSignificanceLevel(ii,jj,1)<0.01 ...
                && mean(psth.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
            test = -1;
        else
            test = 0;
        end
        psth.threeWaysTest(ii,jj,1) = test;
        psth.threeWaysTest_and_salt(ii,jj,1) = test && psth.salt.p_value(ii,jj,1)<0.05;
            
        multipleTest = double([psth.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2)) psth.modulationSignificanceLevel(ii,jj,1)<0.01...
            mean(psth.responsecurveZ(ii,jj,t_duringPulse)) > 1.96 psth.salt.p_value(ii,jj,1)<0.05]);
        multipleTest_sign = -double([psth.rateDuringPulse(ii,jj,1) < ci(2) || isnan(ci(2)) psth.modulationSignificanceLevel(ii,jj,1)<0.01...
            mean(psth.responsecurveZ(ii,jj,t_duringPulse)) < -1.96 psth.salt.p_value(ii,jj,1)<0.05]);
        
        if any(multipleTest_sign([1 3])==-1)
            multipleTest([1 3]) = multipleTest_sign([1 3]);
        end
            
        multipleTest_string = [];
        for zz = 1:length(multipleTest)
            switch multipleTest(zz)
                case -1
                    multipleTest_string = strcat(multipleTest_string, '-');
                case 0
                    multipleTest_string = strcat(multipleTest_string, '=');
                case 1
                    multipleTest_string = strcat(multipleTest_string, '+');
            end
        end
        psth.multipleTest(ii,jj,:) = multipleTest;
        psth.multipleTest_string{ii,jj} = multipleTest_string;

    else
        psth.responsecurve(ii,jj,:) = NaN * psth.responsecurve(ii,jj,:) ;
        psth.responsecurveZ(ii,jj,:) = NaN * psth.responsecurve(ii,jj,:) ;
        psth.modulationSignificanceLevel(ii,jj,1) = NaN;
        psth.rateDuringPulse(ii,jj,1) = NaN;
        psth.rateBeforePulse(ii,jj,1) = NaN;
        psth.rateZDuringPulse(ii,jj,1) = NaN;
         psth.rateZBeforePulse(ii,jj,1) = NaN;
        psth.bootsTrapTest(ii,jj,1) = NaN;
        psth.zscoreTest(ii,jj,1) = NaN;
        psth.threeWaysTest(ii,jj,1) = NaN;
        psth.durationPerPulse(ii,jj,1) = NaN;
        psth.pulseDuration(ii,jj,1) = pulseDuration;
        
        psth.raster.rasterCount{ii,jj} = NaN;
        psth.raster.rasterProb{ii,jj} = NaN;
        psth.raster.TrialsNumber{ii,jj} = NaN;
        psth.raster.times{ii,jj} = NaN;
        psth.raster.rasterTrials{ii,jj} = NaN;
        psth.raster.rasterSpikesTimes{ii,jj} = NaN;
        
        psth.salt.p_value(ii,jj,1) = NaN;
        psth.salt.I_statistics(ii,jj,1) = NaN;
        
        psth.threeWaysTest_and_salt(ii,jj,1) = NaN;
        
        psth.multipleTest(ii,jj,:) = nan(1,4);
        psth.multipleTest_string{ii,jj} = NaN;
    end
    psth.timestamps = t;
end
psth.restricted_intervals = restrict_ints;

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

raster = psth.raster;
if saveMat
    disp('Saving results...');
    raster = psth.raster;
    psth = rmfield(psth,'raster');
    save([basenameFromBasepath(pwd) '.' eventType save_as '.cellinfo.mat'],'psth','-v7.3');
    save([basenameFromBasepath(pwd) '.' eventType save_raster_as '.cellinfo.mat'],'raster','-v7.3');
end

% PLOTS
% 1. Rasters plot
if getRaster && any(any(~isnan(psth.responsecurve)))
    tt = psth.timestamps;
    % st = timestamps;
    % if length(st) > 5000 % if more than 5000
    %     st = randsample(st, 5000);
    %     st = sort(st);
    % end
    disp('   Plotting spikes raster and psth...');
    % [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);
    figure;
    set(gcf,'Position',[200 -500 2500 1200]);
    for jj = 1:size(spikes.UID,2)
        % fprintf(' **Events from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        % rast_x = []; rast_y = [];
        % for kk = 1:length(st)
        %     temp_rast = spikes.times{jj} - st(kk);
        %     temp_rast = temp_rast(temp_rast>winSizePlot(1) & temp_rast<winSizePlot(2));
        %     rast_x = [rast_x temp_rast'];
        %     rast_y = [rast_y kk*ones(size(temp_rast))'];
        % end
        rast_x = raster.rasterSpikesTimes{jj};
        rast_y = raster.rasterTrials{jj};

        % spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
        resp = psth.responsecurveSmooth(jj,:);
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(rast_x, rast_y,'.','MarkerSize',1,'color',[.6 .6 .6])
        hold on
        kk = length(pul);
        plot(tt(tt>raster_time(1) & tt<raster_time(2)), resp(tt>raster_time(1) & tt<raster_time(2)) * kk/max(resp)/2,'k','LineWidth',2);
        xlim([raster_time(1) raster_time(2)]); ylim([0 kk]);

        title([num2str(jj) ' (' psth.multipleTest_string{jj} ')'],'FontWeight','normal','FontSize',10);
        % if psth.threeWaysTest(jj) == 0
        %     title(num2str(jj),'FontWeight','normal','FontSize',10);
        % elseif psth.threeWaysTest(jj) == -1
        %     title([num2str(jj) '(-)'],'FontWeight','normal','FontSize',10);
        % elseif psth.threeWaysTest(jj) == 1
        %     title([num2str(jj) '(+)'],'FontWeight','normal','FontSize',10);
        % end

        if jj == 1
            title([num2str(jj) ' (' psth.multipleTest_string{jj} '), (Bootstrapping, kstest, zscore, salt)'],'FontWeight','normal','FontSize',10);
            ylabel('Event #');
        elseif jj == size(spikes.UID,2)
            xlabel('Time (s)');
        else
            set(gca,'YTick',[],'XTick',[]);
        end
    end
    if savePlot
        if exist('SummaryFigures','dir') == 7
            mkdir('SummaryFigures');
        end
        saveas(gcf,['SummaryFigures\' erase(save_raster_as,'_') eventType '.png']); 
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
    mkdir('SummaryFigures');
    saveas(gcf,['SummaryFigures\' erase(save_as,'_') eventType,'.png']); 
end

cd(prevPath);
end

