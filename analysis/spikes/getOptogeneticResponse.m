
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
% numRep            For boostraping, default, 500. If 0, no bexpoostraping.
% binSize           In seconds, default, 0.001.
% winSize           In seconds, default, 0..
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
% Includes Stimulus-associated spike latency test (SALT) for positive
% responses

% Parse options
p = inputParser;
addParameter(p,'analogChannelsList',NaN);
addParameter(p,'digitalChannelsList',NaN);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',500,@isnumeric);
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
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'salt_time',[-0.250 0.250],@isscalar);
addParameter(p,'salt_win',[0.01],@isscalar);
addParameter(p,'salt_binSize',[0.001],@isscalar);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'onset',0,@isnumeric);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'getRaster',true,@islogical);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','optogeneticResponse',@ischar);
addParameter(p,'save_pulses_as','optogeneticPulses',@ischar);
addParameter(p,'maxNumberOfPulses',5000,@isnumeric);

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
salt_baseline = p.Results.salt_baseline;
salt_time = p.Results.salt_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
bootsTrapCI = p.Results.bootsTrapCI;
onset = p.Results.onset;
offset = p.Results.offset;
getRaster = p.Results.getRaster;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;
save_pulses_as = p.Results.save_pulses_as;
maxNumberOfPulses = p.Results.maxNumberOfPulses;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.optogeneticResponse.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Optogenetic responses already computed! Loading file...');
    load(targetFile.name);
    return
end

ints = [];
session = loadSession;
if isfield(session.epochs{1},'behavioralParadigm') && restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = 'optogeneticResponse_post';
            save_pulses_as = 'optogeneticPulses_post';
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif isfield(session.epochs{1},'behavioralParadigm') && restrict_to_baseline
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

pulsesAnalog.timestamps = []; pulsesAnalog.analogChannelsList = [];

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
if isempty(pulses.timestamps)
    optogeneticResponses = [];
    return 
end

pulses.channel = [pulsesAnalog.analogChannelsList; pulsesDigital.digitalChannelsList + lastAnalogChannels];  % combine pulses
pulses.analogChannelsList = [pulsesAnalog.analogChannelsList; nan(size(pulsesDigital.digitalChannelsList))];  % 
pulses.digitalChannelsList = [nan(size(pulsesAnalog.analogChannelsList)); pulsesDigital.digitalChannelsList];  % 
pulses.duration = round(pulses.timestamps(:,2) - pulses.timestamps(:,1),3);  % 
pulses.isAnalog = [ones(size(pulsesAnalog.analogChannelsList)); zeros(size(pulsesDigital.digitalChannelsList))];
pulses.isDigital = [zeros(size(pulsesAnalog.analogChannelsList)); ones(size(pulsesDigital.digitalChannelsList))];

% restrict_pulses
status = InIntervals(pulses.timestamps(:,1),restrict_ints);
pulses.timestamps = pulses.timestamps(status,:);
pulses.channel = pulses.channel(status,:);
pulses.analogChannelsList = pulses.analogChannelsList(status,:);
pulses.digitalChannelsList = pulses.digitalChannelsList(status,:);
pulses.duration = pulses.duration(status,:);
pulses.isAnalog = pulses.isAnalog(status,:);
pulses.isDigital = pulses.isDigital(status,:);
pulses.restricted_intervals = restrict_ints;

% get cell response
optogeneticResponses = [];
if isempty(pulses.timestamps)
    warning('No pulses were found!');
    return
end

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
pulses.analogChannelsList(toRemove) = [];
pulses.digitalChannelsList(toRemove) = [];
pulses.duration(toRemove) = [];
pulses.isAnalog(toRemove) = [];
pulses.isDigital(toRemove) = [];

%%
spikes = loadSpikes;
% generate random events for boostraping
disp('Computing responses...');

for jj = 1:nConditions
    fprintf('\n Condition %3.i / %3.i \n', jj, nConditions)

    % generate events for boostraping
    disp('Generating boostrap template...');
    nPulses = length(find(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1)));
    randomEvents = [];
    if nPulses > maxNumberOfPulses
        nPulses = maxNumberOfPulses;
    end
    for mm = 1:numRep
        randomEvents{mm} = sort(randsample(timestamps_recording, nPulses))';
    end
    pulseDuration = conditions(jj,1);
    [stccg, t] = CCG([spikes.times randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
    fprintf('\n'); %
    t_duringPulse = t > 0 + onset & t < pulseDuration + offset; 
    randomRatesDuringPulse = squeeze(mean(stccg(t_duringPulse, length(spikes.UID)+1:end,1:length(spikes.UID)),1));
    optogeneticResponses.bootsTrapRate(:,jj) = mean(randomRatesDuringPulse,1);
    optogeneticResponses.bootsTrapRateStd(:,jj) = std(randomRatesDuringPulse,[],1);
    optogeneticResponses.bootsTrapRateSEM(:,jj) = std(randomRatesDuringPulse,[],1)/sqrt(numRep);
    if ~isempty(randomRatesDuringPulse)
        for ii = 1:size(randomRatesDuringPulse,2)
            pd = fitdist(randomRatesDuringPulse(:,ii),'normal');
            optogeneticResponses.bootsTrapCI(ii,jj,1:2) = pd.icdf(bootsTrapCI);
        end
    else
        optogeneticResponses.bootsTrapCI(1:size(randomRatesDuringPulse,2),jj,1:2) = NaN;
    end

    pul = pulses.timestamps(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1),1);
    isAnalog = nanmedian(pulses.isAnalog(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1)));
    channelPulse = median(pulses.channel(pulses.channel == conditions(jj,2) & pulses.duration == conditions(jj,1)));
    if isempty(pulses)
        pul = [0];
    end
    if length(pul)>maxNumberOfPulses
        pul = sort(randsample(pul, maxNumberOfPulses));
    end

    times = spikes.times; times{length(times)+1} = pul;
    [stccg, t] = CCG(times,[],'binSize',binSize,'duration',winSize,'norm','rate'); fprintf('\n'); %
    optogeneticResponses.responsecurve(:,jj,:) = squeeze(stccg(:, end , 1:end-1))';
    if length(times{end}) < minNumberOfPulses
        optogeneticResponses.responsecurve(:,jj,:) = optogeneticResponses.responsecurve(:,jj,:) * NaN;
    end
    t_duringPulse = t > 0 + onset & t < pulseDuration + offset; 
    t_beforePulse = t > -pulseDuration & t < 0; 
    
    numberOfPulses = size(pul,1);
    for ii = 1:size(optogeneticResponses.responsecurve,1)
        if numberOfPulses > minNumberOfPulses
            optogeneticResponses.responsecurveSmooth(ii,jj,:) = smooth(optogeneticResponses.responsecurve(ii,jj,:));
            optogeneticResponses.responsecurveZ(ii,jj,:) = (optogeneticResponses.responsecurve(ii,jj,:)...
                    - mean(optogeneticResponses.responsecurve(ii,jj,t < 0)))...
                    /std(optogeneticResponses.responsecurve(ii,jj,t < 0));
            optogeneticResponses.responsecurveZSmooth(ii,jj,:) = smooth(optogeneticResponses.responsecurveZ(ii,jj,:));
            optogeneticResponses.rateDuringPulse(ii,jj,1) = mean(optogeneticResponses.responsecurve(ii,jj,t_duringPulse));
            optogeneticResponses.rateBeforePulse(ii,jj,1) = mean(optogeneticResponses.responsecurve(ii,jj,t_beforePulse));
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse));
            optogeneticResponses.rateZBeforePulse(ii,jj,1) = mean(optogeneticResponses.responsecurveZ(ii,jj,t_beforePulse));
            optogeneticResponses.durationPerPulse(ii,jj,1) = t(find(t_duringPulse,1,'last')+1) - t(find(t_duringPulse,1,'first')-1);
            optogeneticResponses.pulseDuration(ii,jj,1) = pulseDuration;
            optogeneticResponses.isAnalog(ii,jj,1) = isAnalog;
            optogeneticResponses.isDigital(ii,jj,1) = ~isAnalog;
            optogeneticResponses.channelPulse(ii,jj,1) = channelPulse;
            optogeneticResponses.condition(ii,jj,1) = jj;
            
            try
                [h, optogeneticResponses.modulationSignificanceLevel(ii,jj,1)] = ...
                     kstest2(squeeze(optogeneticResponses.responsecurve(ii,jj,t_duringPulse))...
                        ,squeeze(optogeneticResponses.responsecurve(ii,jj,t_beforePulse)));
            catch
                 optogeneticResponses.modulationSignificanceLevel(ii,jj,1) = NaN;
            end

            % Boostrap test
            ci = squeeze(optogeneticResponses.bootsTrapCI(ii,jj,:));
            if optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2)
                test = 1;
            elseif optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(1)
                test = -1;
            elseif isnan(optogeneticResponses.rateDuringPulse(ii,jj,1))
                    test = NaN;
            else
                test = 0;
            end
            optogeneticResponses.bootsTrapTest(ii,jj,1) = test;

            % z-score change test
            if mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96
                test = 1;
            elseif mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96
                 test = -1;
            elseif isnan(optogeneticResponses.rateZDuringPulse(ii,jj,1))
                test = NaN;
            else
                test = 0;
            end
            optogeneticResponses.zscoreTest(ii,jj,1) = test;
             
            % Generating raster
            rasterX = [];
            rasterY = [];
            if getRaster
                for zz = 1:size(pul,1)
                    temp_spk = spikes.times{ii}(find(spikes.times{ii} - pul(zz,1)  > salt_time(1) & spikes.times{ii} - pul(zz,1)  < salt_time(2))) - pul(zz,1);
                    rasterX = [rasterX; temp_spk];
                    if ~isempty(temp_spk)
                        rasterY = [rasterY; zz * ones(size((temp_spk)))];
                    end
                end
            end
            if ~isempty(rasterX)
                [rasterHist3,c] = hist3([rasterY rasterX],{1:size(pul,1) salt_time(1):salt_binSize:salt_time(2)});
                optogeneticResponses.raster.rasterCount{ii,jj} = rasterHist3;
                optogeneticResponses.raster.rasterProb{ii,jj} = rasterHist3/sum(rasterHist3(:));
                optogeneticResponses.raster.TrialsNumber{ii,jj} = c{1};
                optogeneticResponses.raster.times{ii,jj} = c{2};
                optogeneticResponses.raster.rasterTrials{ii,jj} = rasterY;
                optogeneticResponses.raster.rasterSpikesTimes{ii,jj} = rasterX;

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
                try 
                    [optogeneticResponses.salt.p_value(ii,jj,1), optogeneticResponses.salt.I_statistics(ii,jj,1)] = salt(rasterHist3(:,baseidx(1):baseidx(2)),rasterHist3(:,tidx(1):tidx(2)),salt_binSize, salt_win);
                catch
                    optogeneticResponses.salt.p_value(ii,jj,1) = NaN; 
                    optogeneticResponses.salt.I_statistics(ii,jj,1) = NaN;
                end
            else
                optogeneticResponses.raster.rasterCount{ii,jj} = NaN;
                optogeneticResponses.raster.rasterProb{ii,jj} = NaN;
                optogeneticResponses.raster.TrialsNumber{ii,jj} = NaN;
                optogeneticResponses.raster.times{ii,jj} = NaN;
                optogeneticResponses.raster.rasterTrials{ii,jj} = NaN;
                optogeneticResponses.raster.rasterSpikesTimes{ii,jj} = NaN;
                optogeneticResponses.salt.p_value(ii,jj,1) = NaN;
                optogeneticResponses.salt.I_statistics(ii,jj,1) = NaN;
            end

            % multiple test test. If not boostrap, it would be 2 ways.
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
            optogeneticResponses.threeWaysTest_and_salt(ii,jj,1) = test && optogeneticResponses.salt.p_value(ii,jj,1)<0.05;
            
            multipleTest = double([optogeneticResponses.rateDuringPulse(ii,jj,1) > ci(2) || isnan(ci(2)) optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01...
                mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) > 1.96 optogeneticResponses.salt.p_value(ii,jj,1)<0.05]);
            multipleTest_sign = -double([optogeneticResponses.rateDuringPulse(ii,jj,1) < ci(2) || isnan(ci(2)) optogeneticResponses.modulationSignificanceLevel(ii,jj,1)<0.01...
                mean(optogeneticResponses.responsecurveZ(ii,jj,t_duringPulse)) < -1.96 optogeneticResponses.salt.p_value(ii,jj,1)<0.05]);
            
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
            optogeneticResponses.multipleTest(ii,jj,:) = multipleTest;
            optogeneticResponses.multipleTest_string{ii,jj} = multipleTest_string;
        else
            optogeneticResponses.responsecurve(ii,jj,:) = NaN * optogeneticResponses.responsecurve(ii,jj,:) ;
            optogeneticResponses.responsecurveZ(ii,jj,:) = NaN * optogeneticResponses.responsecurve(ii,jj,:) ;
            optogeneticResponses.modulationSignificanceLevel(ii,jj,1) = NaN;
            optogeneticResponses.rateDuringPulse(ii,jj,1) = NaN;
            optogeneticResponses.rateBeforePulse(ii,jj,1) = NaN;
            optogeneticResponses.rateZDuringPulse(ii,jj,1) = NaN;
             optogeneticResponses.rateZBeforePulse(ii,jj,1) = NaN;
            optogeneticResponses.bootsTrapTest(ii,jj,1) = NaN;
            optogeneticResponses.zscoreTest(ii,jj,1) = NaN;
            optogeneticResponses.threeWaysTest(ii,jj,1) = NaN;
            optogeneticResponses.durationPerPulse(ii,jj,1) = NaN;
            optogeneticResponses.isAnalog(ii,jj,1) = isAnalog;
            optogeneticResponses.isDigital(ii,jj,1) = isAnalog - 1;
            optogeneticResponses.channelPulse(ii,jj,1) = channelPulse;
            optogeneticResponses.pulseDuration(ii,jj,1) = pulseDuration;
            
            optogeneticResponses.raster.rasterCount{ii,jj} = NaN;
            optogeneticResponses.raster.rasterProb{ii,jj} = NaN;
            optogeneticResponses.raster.TrialsNumber{ii,jj} = NaN;
            optogeneticResponses.raster.times{ii,jj} = NaN;
            optogeneticResponses.raster.rasterTrials{ii,jj} = NaN;
            optogeneticResponses.raster.rasterSpikesTimes{ii,jj} = NaN;
            
            optogeneticResponses.salt.p_value(ii,jj,1) = NaN;
            optogeneticResponses.salt.I_statistics(ii,jj,1) = NaN;
            
            optogeneticResponses.threeWaysTest_and_salt(ii,jj,1) = NaN;
            
            optogeneticResponses.multipleTest(ii,jj,:) = nan(1,4);
            optogeneticResponses.multipleTest_string{ii,jj} = NaN;
        end
        optogeneticResponses.timestamps = t;   
    end
end
optogeneticResponses.restricted_intervals = restrict_ints;

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

optogeneticResponses.salt_baseline = salt_baseline;
optogeneticResponses.salt_time = salt_time;
optogeneticResponses.salt_win = salt_win;
optogeneticResponses.salt_binSize = salt_binSize;

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
    optogeneticResponses_raster = optogeneticResponses;
    optogeneticResponses = rmfield(optogeneticResponses,'raster');
    disp(' Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.' save_as '.cellinfo.mat'],'optogeneticResponses','-v7.3');
    if getRaster
        save([basenameFromBasepath(pwd) '.' save_as '_raster.cellinfo.mat'],'optogeneticResponses_raster','-v7.3');
    end
end

if saveEventsFile
    disp('Saving timestamps...');
    filename = split(pwd,filesep); filename = filename{end};
    optoPulses = optogeneticResponses.pulses;
    optoPulses.stimulationEpochs = stimulationEpochs;
    save([filename '.' save_pulses_as '.events.mat'],'optoPulses');
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
            
            title([num2str(jj) ' (' optogeneticResponses.multipleTest_string{jj,ii} ')'],'FontWeight','normal','FontSize',10);

            if jj == 1
                title([num2str(jj) ' (' optogeneticResponses.multipleTest_string{jj,ii} '), (Bootstrapping, kstest, zscore, salt)'],'FontWeight','normal','FontSize',10);
                ylabel('Trial');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\' save_as 'Raster_ch',num2str(conditions(ii,2)),'_dur',num2str(conditions(ii,1)),'ch.png']); 
%         close(gcf);
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
    saveas(gcf,['SummaryFigures\' save_as 'Psth.png']); 
end          


cd(prevPath);
end