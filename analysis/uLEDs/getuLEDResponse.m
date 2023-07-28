
function [uLEDResponses] = getuLEDResponse(varargin)
% [uLEDResponses] = getuLEDResponse(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% during uLED stimulation
%
% <OPTIONALS>
% uLEDPulses        uLEDPulses structure, output from getuLEDPulses.
% spikes            buzcode spikes structure, if not provided tries loadSpikes.
% basepath          By default pwd.
% numRep            For boostraping, default, 500. If 0, no boostraping.
% binSize           In seconds, default, 0.001.
% winSize           In seconds, default, 0.5.
% offset            Numeric modifier for the end of the time window
%                       that will be use for assesing neurons responses,
%                       default [0].Use array to specifiy diferent onsets
%                       for different pulse conditions (Ex. [10 0])       
% onset             Numeric modifier for the beggining of the time window
%                       that will be use for assesing neurons responses,
%                       default [0]. Use array to specifiy diferent onsets
%                       for different pulse conditions.
% doPlot            Default true.
% winSizePlot       Default [-0.1 .5];
% force             Default, false.                   
%
% OUTPUTS
% uLEDResponses
%
% Manu Valero 2022
% Includes Stimulus-associated spike latency test (SALT) for positive
% responses

% Parse options
p = inputParser;
addParameter(p,'uLEDPulses',NaN);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',500,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',.1,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'onset',0,@isnumeric);
addParameter(p,'winSizePlot',[-.02 .05],@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'minNumberOfPulses',200,@isnumeric);
addParameter(p,'duration_round_decimal',3,@isscalar);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'salt_time',[-0.250 0.250],@isscalar);
addParameter(p,'salt_win',[0.005],@isscalar);
addParameter(p,'salt_binSize',[0.001],@isscalar);

parse(p, varargin{:});
uLEDPulses = p.Results.uLEDPulses;
basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
doPlot = p.Results.doPlot;
offset = p.Results.offset;
onset = p.Results.onset;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;
minNumberOfPulses = p.Results.minNumberOfPulses;
duration_round_decimal = p.Results.duration_round_decimal;
salt_baseline = p.Results.salt_baseline;
salt_time = p.Results.salt_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
bootsTrapCI = p.Results.bootsTrapCI;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDResponse.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('uLED responses already computed! Loading file...');
    load(targetFile.name);
    return
end

if isnan(uLEDPulses)
    uLEDPulses = getuLEDPulses;
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

% Get cell responses!!
uLEDResponses = [];
codes = 1:max(uLEDPulses.code);
timestamps_recording = min(uLEDPulses.timestamps(:,2)):1/1250:max(uLEDPulses.timestamps(:,2));

if length(uLEDPulses.conditionDurationID)>length(onset) && length(onset)==1
    onset(1:length(uLEDPulses.conditionDurationID)) = onset(1);
end
if length(uLEDPulses.conditionDurationID)>length(offset) && length(offset)==1
    offset(1:length(uLEDPulses.conditionDurationID)) = offset(1);
end

for kk = 1:length(uLEDPulses.conditionDurationID)
    fprintf('> Condition %3.i/ %3.i \n',kk, length(uLEDPulses.conditionDurationID)); %
    
    pulseDuration = uLEDPulses.conditionDuration(kk);
    % generate random events for boostraping
    nPulses = int32(length(find(uLEDPulses.conditionID == uLEDPulses.conditionDurationID(kk)))/...
            length(unique(uLEDPulses.code(uLEDPulses.conditionID == uLEDPulses.conditionDurationID(kk)))));
    randomEvents = [];
    disp('Generating boostrap template...');
    for mm = 1:numRep
        randomEvents{mm} = sort(randsample(timestamps_recording, nPulses))';
    end
    % loop for spike response
    disp('Computing cell responses... ');
    for ii = 1:length(spikes.UID)
        fprintf(' **Pulses from unit %3.i/ %3.i \n',ii, size(spikes.UID,2)); %
        [stccg, t] = CCG([spikes.times{ii} randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
        t_duringPulse = t > 0 + onset(kk) & t < pulseDuration + offset(kk); 
        randomRatesDuringPulse = mean(stccg(t_duringPulse,2:size(randomEvents,2)+1,1),1);
        uLEDResponses.bootsTrapRate(ii,kk) = mean(randomRatesDuringPulse);
        uLEDResponses.bootsTrapRateStd(ii,kk) = std(randomRatesDuringPulse);
        pd = fitdist(randomRatesDuringPulse','normal');
        uLEDResponses.bootsTrapCI(ii,kk,:) = pd.icdf(bootsTrapCI);
        for jj = 1:length(codes)
            pulses = uLEDPulses.timestamps(uLEDPulses.code == codes(jj) & uLEDPulses.conditionID==kk,:);
            if length(pulses) > minNumberOfPulses
                [stccg, t] = CCG({spikes.times{ii}, pulses(:,1)},[],'binSize',binSize,'duration',winSize,'norm','rate');
                uLEDResponses.responsecurve(ii,kk,jj,:) = stccg(:,2,1);
                uLEDResponses.responsecurveSmooth(ii,kk,jj,:) = smooth(stccg(:,2,1));
                t_duringPulse = t > 0 + onset(kk) & t < pulseDuration + offset(kk); 
                t_beforePulse = t > -pulseDuration & t < 0; 

                uLEDResponses.responsecurveZ(ii,kk,jj,:) = (stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1));
                uLEDResponses.responsecurveZSmooth(ii,kk,jj,:) = smooth(stccg(:,2,1) - mean(stccg(t < 0,2,1)))/std(stccg(t < 0,2,1));
                uLEDResponses.rateDuringPulse(ii,kk,jj,1) = mean(stccg(t_duringPulse,2,1));
                uLEDResponses.rateBeforePulse(ii,kk,jj,1) = mean(stccg(t_beforePulse,2,1));
                uLEDResponses.rateZDuringPulse(ii,kk,jj,1) = mean(squeeze(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)));
                uLEDResponses.rateZBeforePulse(ii,kk,jj,1) = mean(squeeze(uLEDResponses.responsecurveZ(ii,kk,jj,t_beforePulse)));
                uLEDResponses.codes(ii,kk,jj,1) = codes(jj);
                [h, uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1)] = kstest2(stccg(t_duringPulse,2,1),stccg(t_beforePulse,2,1));

                % Boostrap test
                ci = squeeze(uLEDResponses.bootsTrapCI(ii,kk,:));
                if uLEDResponses.rateDuringPulse(ii,kk,jj,1) > ci(2)
                    test = 1;
                elseif uLEDResponses.rateDuringPulse(ii,kk,jj,1) < ci(1)
                    test = -1;
                else
                    test = 0;
                end
                uLEDResponses.bootsTrapTest(ii,kk,jj,1) = test;

                % z-score change test
                if mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) > 1.96
                    test = 1;
                elseif mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) < -1.96
                    test = -1;
                else
                    test = 0;
                end
                uLEDResponses.zscoreTest(ii,kk,jj,1) = test;

                % Generating raster
                rasterX = [];
                rasterY = [];
                for zz = 1:size(pulses,1)
                    temp_spk = spikes.times{ii}(find(spikes.times{ii} - pulses(zz,1)  > salt_time(1) & spikes.times{ii} - pulses(zz,1)  < salt_time(2))) - pulses(zz,1);
                    rasterX = [rasterX; temp_spk];
                    if ~isempty(temp_spk)
                        rasterY = [rasterY; zz * ones(size((temp_spk)))];
                    end
                end
                if ~isempty(rasterX)
                    [rasterHist3,c] = hist3([rasterY rasterX],{1:size(pulses,1) salt_time(1):salt_binSize:salt_time(2)});
                    uLEDResponses.raster.rasterCount{ii,kk,jj} = rasterHist3;
                    uLEDResponses.raster.rasterProb{ii,kk,jj} = rasterHist3/sum(rasterHist3(:));
                    uLEDResponses.raster.TrialsNumber{ii,kk,jj} = c{1};
                    uLEDResponses.raster.times{ii,kk,jj} = c{2};
                    uLEDResponses.raster.rasterTrials{ii,kk,jj} = rasterY;
                    uLEDResponses.raster.rasterSpikesTimes{ii,kk,jj} = rasterX;

                    time = c{2};
                    baseidx = dsearchn(time', salt_baseline');
                    tidx = dsearchn(time', [0; pulseDuration*2]);    
                    st = length(baseidx(1):baseidx(2));
                    nmbn = round(salt_win/salt_binSize);
                    v = 1:nmbn:st;
                    if any((v + nmbn - 1) > st)
                        error('reduce window size or baseline duration')
                    end
                    [uLEDResponses.salt.p_value(ii,kk,jj,1), uLEDResponses.salt.I_statistics(ii,kk,jj,1)] = salt(rasterHist3(:,baseidx(1):baseidx(2)),rasterHist3(:,tidx(1):tidx(2)),salt_binSize, salt_win);
                else
                    uLEDResponses.raster.rasterCount{ii,kk,jj} = NaN;
                    uLEDResponses.raster.rasterProb{ii,kk,jj} = NaN;
                    uLEDResponses.raster.TrialsNumber{ii,kk,jj} = NaN;
                    uLEDResponses.raster.times{ii,kk,jj} = NaN;
                    uLEDResponses.raster.rasterTrials{ii,kk,jj} = NaN;
                    uLEDResponses.raster.rasterSpikesTimes{ii,kk,jj} = NaN;
                    uLEDResponses.salt.p_value(ii,kk,jj,1) = NaN;
                    uLEDResponses.salt.I_statistics(ii,kk,jj,1) = NaN;
                end
                % multiple test test. If not boostrap, it would be 2 ways.
                if (uLEDResponses.rateDuringPulse(ii,kk,jj,1) > ci(2) || isnan(ci(2))) && uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1)<0.01...
                        && mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) > 1.96
                    test = 1;
                elseif (uLEDResponses.rateDuringPulse(ii,kk,jj,1) < ci(1) || isnan(ci(1))) && uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1)<0.01 ...
                        && mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) < -1.96
                    test = -1;
                else
                    test = 0;
                end
                uLEDResponses.threeWaysTest(ii,kk,jj,1) = test;
                uLEDResponses.threeWaysTest_and_salt(ii,kk,jj,1) = int32(test && uLEDResponses.salt.p_value(ii,kk,jj,1)<0.05);

                multipleTest = double([uLEDResponses.rateDuringPulse(ii,kk,jj,1) > ci(2) || isnan(ci(2)) uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1)<0.01...
                    mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) > 1.96 uLEDResponses.salt.p_value(ii,kk,jj,1)<0.05]);
                multipleTest_sign = -double([uLEDResponses.rateDuringPulse(ii,kk,jj,1) < ci(2) || isnan(ci(2)) uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1)<0.01...
                    mean(uLEDResponses.responsecurveZ(ii,kk,jj,t_duringPulse)) < -1.96 uLEDResponses.salt.p_value(ii,kk,jj,1)<0.05]);

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
                uLEDResponses.multipleTest(ii,kk,jj,:) = multipleTest;
                uLEDResponses.multipleTest_string{ii,kk,jj} = multipleTest_string;
                
                uLEDResponses.numberOfPulses(ii,kk,jj,1) = size(pulses,1);
                uLEDResponses.pulseDuration(ii,kk,jj,1) = uLEDPulses.conditionDuration(kk);
            else
                uLEDResponses.responsecurve(ii,kk,jj,:) = nan(winSize/binSize + 1,1);
                uLEDResponses.responsecurveSmooth(ii,kk,jj,:) = nan(winSize/binSize + 1,1);

                uLEDResponses.responsecurveZ(ii,kk,jj,:) = nan(winSize/binSize + 1,1);
                uLEDResponses.responsecurveZSmooth(ii,kk,jj,:) = nan(winSize/binSize + 1,1);
                uLEDResponses.rateDuringPulse(ii,kk,jj,1) = NaN;
                uLEDResponses.rateBeforePulse(ii,kk,jj,1) = NaN;
                uLEDResponses.rateZDuringPulse(ii,kk,jj,1) = NaN;
                uLEDResponses.rateZBeforePulse(ii,kk,jj,1) = NaN;
                uLEDResponses.codes(ii,kk,jj,1) =NaN;
                uLEDResponses.modulationSignificanceLevel(ii,kk,jj,1) = NaN;

                uLEDResponses.bootsTrapTest(ii,kk,jj,1) = NaN;
                uLEDResponses.zscoreTest(ii,kk,jj,1) = NaN;

                uLEDResponses.raster.rasterCount{ii,kk,jj} = NaN;
                uLEDResponses.raster.rasterProb{ii,kk,jj} = NaN;
                uLEDResponses.raster.TrialsNumber{ii,kk,jj} = NaN;
                uLEDResponses.raster.times{ii,kk,jj} = NaN;
                uLEDResponses.raster.rasterTrials{ii,kk,jj} = NaN;
                uLEDResponses.raster.rasterSpikesTimes{ii,kk,jj} = NaN;
                uLEDResponses.salt.p_value(ii,kk,jj,1) = NaN;
                uLEDResponses.salt.I_statistics(ii,kk,jj,1) = NaN;

                uLEDResponses.threeWaysTest(ii,kk,jj,1) = NaN;
                uLEDResponses.threeWaysTest_and_salt(ii,kk,jj,1) = NaN;

                uLEDResponses.multipleTest(ii,kk,jj,:) = NaN;
                uLEDResponses.multipleTest_string{ii,kk,jj} = NaN;
                
                uLEDResponses.numberOfPulses(ii,kk,jj,1) = NaN;
                uLEDResponses.pulseDuration(ii,kk,jj,1) =NaN;
            end
        end
        uLEDResponses.timestamps = t;
    end
end

% parse cell responses
for kk = 1:length(uLEDPulses.conditionDurationID)
    for ii = 1:length(spikes.UID)
        uLEDResponses.noRespLEDs.LEDs{ii,kk} = find(uLEDResponses.bootsTrapTest(ii,kk,:) == 0);
        uLEDResponses.respLEDs.LEDs{ii,kk} = find(uLEDResponses.bootsTrapTest(ii,kk,:) == 1);
        uLEDResponses.negRespLEDs.LEDs{ii,kk} = find(uLEDResponses.bootsTrapTest(ii,kk,:)== -1);

        ratio = squeeze(uLEDResponses.rateDuringPulse(ii,kk,:)./uLEDResponses.rateBeforePulse(ii,kk,:));
        ratio(isinf(ratio)) = NaN;
        maxRespLED = find(max(ratio)==ratio,1);
        minRespLED = find(min(ratio)==ratio,1);         
        if uLEDResponses.bootsTrapTest(ii,kk,maxRespLED) == 1
            uLEDResponses.maxRespLED.LEDs(ii,kk) = maxRespLED;
        else
            uLEDResponses.maxRespLED.LEDs(ii,kk) = NaN;
        end

        if uLEDResponses.bootsTrapTest(ii,kk,minRespLED) == -1
            uLEDResponses.minRespLED.LEDs(ii,kk) = minRespLED;
        else
            uLEDResponses.minRespLED.LEDs(ii,kk) = NaN;
        end

        % noRespLEDs
        uLEDResponses.noRespLEDs.values{ii,kk} = uLEDResponses.rateDuringPulse(ii,kk,uLEDResponses.noRespLEDs.LEDs{ii,kk});
        uLEDResponses.noRespLEDs.valuesBeforePulse{ii,kk} = uLEDResponses.rateBeforePulse(ii,kk,uLEDResponses.noRespLEDs.LEDs{ii,kk});
        uLEDResponses.noRespLEDs.ratioBeforeAfter{ii,kk} = uLEDResponses.noRespLEDs.values{ii,kk}./uLEDResponses.noRespLEDs.valuesBeforePulse{ii,kk}; 
        uLEDResponses.noRespLEDs.meanRatio(ii,kk) = mean(uLEDResponses.noRespLEDs.ratioBeforeAfter{ii,kk});
        uLEDResponses.noRespLEDs.meanRateBeforePulse(ii,kk) = mean(uLEDResponses.noRespLEDs.valuesBeforePulse{ii,kk});
        uLEDResponses.noRespLEDs.meanRateDuringPulse(ii,kk) = mean(uLEDResponses.noRespLEDs.values{ii,kk});
        if isnan(uLEDResponses.noRespLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
            uLEDResponses.noRespLEDs.LEDs{ii,kk} = [];
            uLEDResponses.noRespLEDs.values{ii,kk} = [];
            uLEDResponses.noRespLEDs.valuesBeforePulse{ii,kk} = [];
            uLEDResponses.noRespLEDs.ratioBeforeAfter{ii,kk} = [];
        end

        % negRespLEDs
        uLEDResponses.negRespLEDs.values{ii,kk} = uLEDResponses.rateDuringPulse(ii,kk,uLEDResponses.negRespLEDs.LEDs{ii,kk});
        uLEDResponses.negRespLEDs.valuesBeforePulse{ii,kk} = uLEDResponses.rateBeforePulse(ii,kk,uLEDResponses.negRespLEDs.LEDs{ii,kk});  
        uLEDResponses.negRespLEDs.ratioBeforeAfter{ii,kk} = uLEDResponses.negRespLEDs.values{ii,kk}./uLEDResponses.negRespLEDs.valuesBeforePulse{ii,kk};  
        uLEDResponses.negRespLEDs.meanRatio(ii,kk) = mean(uLEDResponses.negRespLEDs.ratioBeforeAfter{ii,kk});
        uLEDResponses.negRespLEDs.ratioNoResp(ii,kk) = mean(uLEDResponses.negRespLEDs.values{ii,kk})./nanmean(uLEDResponses.noRespLEDs.values{ii,kk});
        uLEDResponses.negRespLEDs.meanRateBeforePulse(ii,kk) = mean(uLEDResponses.negRespLEDs.valuesBeforePulse{ii,kk});
        uLEDResponses.negRespLEDs.meanRateDuringPulse(ii,kk) = mean(uLEDResponses.negRespLEDs.values{ii,kk});
        if isnan(uLEDResponses.negRespLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
            uLEDResponses.negRespLEDs.LEDs{ii,kk} = [];
            uLEDResponses.negRespLEDs.values{ii,kk} = [];
            uLEDResponses.negRespLEDs.valuesBeforePulse{ii,kk} = [];
            uLEDResponses.negRespLEDs.ratioBeforeAfter{ii,kk} = [];
        end

        % respLEDs
        uLEDResponses.respLEDs.values{ii,kk} = uLEDResponses.rateDuringPulse(ii,kk,uLEDResponses.respLEDs.LEDs{ii,kk});
        uLEDResponses.respLEDs.valuesBeforePulse{ii,kk} = uLEDResponses.rateBeforePulse(ii,kk,uLEDResponses.respLEDs.LEDs{ii,kk});  
        uLEDResponses.respLEDs.ratioBeforeAfter{ii,kk} = uLEDResponses.respLEDs.values{ii}./uLEDResponses.respLEDs.valuesBeforePulse{ii,kk};  
        uLEDResponses.respLEDs.meanRatio(ii,kk) = mean(uLEDResponses.respLEDs.ratioBeforeAfter{ii,kk});
        uLEDResponses.respLEDs.ratioNoResp(ii,kk) = nanmean(uLEDResponses.respLEDs.values{ii,kk})./nanmean(uLEDResponses.noRespLEDs.values{ii,kk});
        uLEDResponses.respLEDs.meanRateBeforePulse(ii,kk) = mean(uLEDResponses.respLEDs.valuesBeforePulse{ii,kk});
        uLEDResponses.respLEDs.meanRateDuringPulse(ii,kk) = mean(uLEDResponses.respLEDs.values{ii,kk});
        if isnan(uLEDResponses.respLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
            uLEDResponses.respLEDs.LEDs{ii,kk} = [];
            uLEDResponses.respLEDs.values{ii,kk} = [];
            uLEDResponses.respLEDs.valuesBeforePulse{ii,kk} = [];
            uLEDResponses.respLEDs.ratioBeforeAfter{ii,kk} = [];
        end

        % maxRespLED
        if ~isnan(uLEDResponses.maxRespLED.LEDs(ii,kk))
            uLEDResponses.maxRespLED.values(ii,kk) = uLEDResponses.rateDuringPulse(ii,kk,uLEDResponses.maxRespLED.LEDs(ii,kk));
            uLEDResponses.maxRespLED.valuesBeforePulse(ii,kk) = uLEDResponses.rateBeforePulse(ii,kk,uLEDResponses.maxRespLED.LEDs(ii,kk));   
            uLEDResponses.maxRespLED.ratioBeforeAfter(ii,kk) = uLEDResponses.maxRespLED.values(ii,kk)./uLEDResponses.maxRespLED.valuesBeforePulse(ii,kk); 
            uLEDResponses.maxRespLED.ratioNoResp(ii,kk) = uLEDResponses.maxRespLED.values(ii)./nanmean(uLEDResponses.noRespLEDs.values{ii,kk});
            if isnan(uLEDResponses.maxRespLED.values(ii,kk))
                uLEDResponses.maxRespLED.LEDs(ii,kk) = NaN;
            end
        else
            uLEDResponses.maxRespLED.values(ii,kk) = NaN;
            uLEDResponses.maxRespLED.valuesBeforePulse(ii,kk) = NaN;
            uLEDResponses.maxRespLED.ratioBeforeAfter(ii,kk) = NaN;
            uLEDResponses.maxRespLED.ratioNoResp(ii,kk) = NaN;
        end
        
        % minRespLED
        if ~isnan(uLEDResponses.minRespLED.LEDs(ii,kk))
            uLEDResponses.minRespLED.values(ii,kk) = uLEDResponses.rateDuringPulse(ii,kk,uLEDResponses.minRespLED.LEDs(ii,kk));
            uLEDResponses.minRespLED.valuesBeforePulse(ii,kk) = uLEDResponses.rateBeforePulse(ii,kk,uLEDResponses.minRespLED.LEDs(ii,kk));   
            uLEDResponses.minRespLED.ratioBeforeAfter(ii,kk) = uLEDResponses.minRespLED.values(ii,kk)./uLEDResponses.minRespLED.valuesBeforePulse(ii,kk);  
            uLEDResponses.minRespLED.ratioNoResp(ii,kk) = uLEDResponses.minRespLED.values(ii,kk)./nanmean(uLEDResponses.noRespLEDs.values{ii,kk});
            if isnan(uLEDResponses.minRespLED.values(ii,kk))
                uLEDResponses.minRespLED.LEDs(ii,kk) = NaN;
            end
        else
            uLEDResponses.minRespLED.values(ii,kk) = NaN;
            uLEDResponses.minRespLED.valuesBeforePulse(ii,kk) = NaN;   
            uLEDResponses.minRespLED.ratioBeforeAfter(ii,kk) = NaN;  
            uLEDResponses.minRespLED.ratioNoResp(ii,kk) = NaN;
        end

    end
end
uLEDResponses.ratioBeforeAfter = uLEDResponses.rateDuringPulse./uLEDResponses.rateBeforePulse;

% find cells in non-stimulated shanks
uLEDResponses.unisInNonStimulatedShanks = (spikes.shankID == uLEDPulses.nonStimulatedShank);

% parse non-responsive and responsive cells, only in group 1 (no manipulated)
uLEDResponses.drivenCells = [];
for kk = 1:length(uLEDPulses.conditionDurationID)
    for ii = 1:length(spikes.UID)
        if any(uLEDResponses.bootsTrapTest(ii,kk,:)==1)
            uLEDResponses.drivenCells(ii,kk) = 1;
        else
            uLEDResponses.drivenCells(ii,kk) = 0;
        end

        if any(uLEDResponses.bootsTrapTest(ii,kk,:)==-1)
            uLEDResponses.inhibitedCells(ii,kk) = 1;
        else
            uLEDResponses.inhibitedCells(ii,kk) = 0;
        end
    end
end

uLEDResponses_raster = uLEDResponses;
uLEDResponses = rmfield(uLEDResponses,'raster');
disp('Saving...');
save([basenameFromBasepath(pwd) '.uLEDResponse.cellinfo.mat'],'uLEDResponses');
save([basenameFromBasepath(pwd) '.uLEDResponse_raster.cellinfo.mat'],'uLEDResponses_raster','-v7.3');

if doPlot
    statisticDots = linspace(-0.015, -0.005, 4);
    nLEDS = size(uLEDResponses.responsecurve,3);
    for kk = 1:length(uLEDPulses.conditionDurationID)
        figure;
        set(gcf,'Position',[100 -600 2500 1200]);
        tiledlayout(10,ceil(size(spikes.UID,2)/10),'TileSpacing','tight','Padding','tight');
        for ii = 1:size(uLEDResponses.bootsTrapCI,1)
            % subplot(10,ceil(size(spikes.UID,2)/10),ii); %
            nexttile
            imagesc(uLEDResponses.timestamps, 1:nLEDS, squeeze(uLEDResponses.responsecurve(ii,kk,:,:))); caxis([0 30]); 
            xlim(winSizePlot); ylim([-0.5 nLEDS+0.5])
            hold on
            plot([0 0], [0 nLEDS+2],'w','LineWidth',1.5);
            plot([uLEDResponses.pulseDuration(ii,kk,1) uLEDResponses.pulseDuration(ii,kk,1)], [0 nLEDS+2],'color',[.5 .5 .5],'LineWidth',1.5);
            title(num2str(ii),'FontWeight','normal','FontSize',10);
            set(gca,'TickDir','out')
            
            if ii == 1
                ylabel('LED [#][s][0-30Hz]');
            end
            
            modulationSignificanceLevel = squeeze(uLEDResponses.modulationSignificanceLevel(ii,kk,:));
            bootsTrapTest = squeeze(uLEDResponses.bootsTrapTest(ii,kk,:));
            for mm = 1:length(modulationSignificanceLevel)
                if bootsTrapTest(mm,1) == 1
                    if modulationSignificanceLevel(mm) < 0.05 && modulationSignificanceLevel(mm) > 0.01
                        plot(statisticDots(1), [mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.01 && modulationSignificanceLevel(mm) > 0.001
                        plot(statisticDots(1:2), [mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.001 && modulationSignificanceLevel(mm) > 0.0001
                        plot(statisticDots(1:3), [mm mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.0001
                        plot(statisticDots(1:4), [mm mm mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    end
                end
            end
        end
        exportgraphics(gcf,['SummaryFigures\uLEDResponse_rate_condition',num2str(kk),'_dur',num2str(uLEDResponses.pulseDuration(ii,kk,1)),'s.png']);
        
        figure;
        set(gcf,'Position',[100 -600 2500 1200]);
        tiledlayout(10,ceil(size(spikes.UID,2)/10),'TileSpacing','tight','Padding','tight');
        for ii = 1:size(uLEDResponses.bootsTrapCI,1)
            % subplot(10,ceil(size(spikes.UID,2)/10),ii); %
            nexttile
            imagesc(uLEDResponses.timestamps, 1:nLEDS, squeeze(uLEDResponses.responsecurveZ(ii,kk,:,:))); caxis([-5 5]); 
            xlim(winSizePlot); ylim([-0.5 nLEDS+0.5])
            hold on
            plot([0 0], [0 nLEDS+2],'w','LineWidth',1.5);
            plot([uLEDResponses.pulseDuration(ii,kk,1) uLEDResponses.pulseDuration(ii,kk,1)], [0 nLEDS+2],'color',[.5 .5 .5],'LineWidth',1.5);
            title(num2str(ii),'FontWeight','normal','FontSize',10);
            set(gca,'TickDir','out')
            
            if ii == 1
                ylabel('LED [#][s][+-5 s.d.]');
            end
            
            modulationSignificanceLevel = squeeze(uLEDResponses.modulationSignificanceLevel(ii,kk,:));
            bootsTrapTest = squeeze(uLEDResponses.bootsTrapTest(ii,kk,:));
            for mm = 1:length(modulationSignificanceLevel)
                if bootsTrapTest(mm,1) == 1
                    if modulationSignificanceLevel(mm) < 0.05 && modulationSignificanceLevel(mm) > 0.01
                        plot(statisticDots(1), [mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.01 && modulationSignificanceLevel(mm) > 0.001
                        plot(statisticDots(1:2), [mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.001 && modulationSignificanceLevel(mm) > 0.0001
                        plot(statisticDots(1:3), [mm mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    elseif modulationSignificanceLevel(mm) < 0.0001
                        plot(statisticDots(1:4), [mm mm mm mm],'ow', 'MarkerFaceColor','w','MarkerSize',3);
                    end
                end
            end
        end
        exportgraphics(gcf,['SummaryFigures\uLEDResponse_Zscored_condition',num2str(kk),'_dur',num2str(uLEDResponses.pulseDuration(ii,kk,1)),'s.png']);
    end
end

cd(prevPath);
end
