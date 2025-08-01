
function [uLEDResponses_interval] = getuLEDResponse_intervals(intervals,varargin)
% [uLEDResponses] = getuLEDResponse_intervals(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% during uLED stimulation that ocurr (or not) at a given interval
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
% Manu Valero 2023

% Parse options
p = inputParser;
addRequired(p,'intervals',@isnumeric);
addParameter(p,'uLEDPulses',NaN);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',500,@isnumeric);
addParameter(p,'numRep_fake',500,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',.1,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'getRaster',true,@islogical);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'onset',0,@isnumeric);
addParameter(p,'winSizePlot',[-.02 .05],@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'minNumberOfPulses',200,@isnumeric);
% addParameter(p,'duration_round_decimal',3,@isscalar);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'salt_time',[-0.250 0.250],@isscalar);
addParameter(p,'salt_win',0.005,@isscalar);
addParameter(p,'salt_binSize',0.001,@isscalar);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'duration_pulse',0.020,@isnumeric);
addParameter(p,'save_as','uLEDResponse_interval',@ischar);
addParameter(p,'boostraping_type','pulses',@ischar);
addParameter(p,'interpolate_pulse_sides',true,@islogical); % 

% declare global variables
global salt_baseline salt_time salt_win salt_binSize

parse(p, intervals, varargin{:});
uLEDPulses = p.Results.uLEDPulses;
basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
numRep_fake = p.Results.numRep_fake;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
doPlot = p.Results.doPlot;
offset = p.Results.offset;
onset = p.Results.onset;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;
minNumberOfPulses = p.Results.minNumberOfPulses;
% duration_round_decimal = p.Results.duration_round_decimal;
salt_baseline = p.Results.salt_baseline;
salt_time = p.Results.salt_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
bootsTrapCI = p.Results.bootsTrapCI;
getRaster = p.Results.getRaster;
verbose = p.Results.verbose;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
duration_pulse = p.Results.duration_pulse;
save_as = p.Results.save_as;
boostraping_type = p.Results.boostraping_type;
interpolate_pulse_sides = p.Results.interpolate_pulse_sides;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDResponse_interval.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('uLED responses already computed! Loading file...');
    load(targetFile.name);
    return
end

if ~isstruct(uLEDPulses) && isnan(uLEDPulses)
    uLEDPulses = getuLEDPulses;
end
% remove unncesary fields
uLEDPulses_temp = uLEDPulses;
clear uLEDPulses
uLEDPulses.conditionID = uLEDPulses_temp.conditionID;
uLEDPulses.list_of_conditions = uLEDPulses_temp.list_of_conditions;
uLEDPulses.list_of_durations = uLEDPulses_temp.list_of_durations;
uLEDPulses.list_of_epochs = uLEDPulses_temp.list_of_epochs;
uLEDPulses.code = uLEDPulses_temp.code;
uLEDPulses.timestamps = uLEDPulses_temp.timestamps;
uLEDPulses.nonStimulatedShank = uLEDPulses_temp.nonStimulatedShank;
clear uLEDPulses_temp

%check that the pulse duration is 0.02 
CONDITION =[];

if length(uLEDPulses.list_of_conditions) > 1
       
    CONDITION_LIST = uLEDPulses.list_of_conditions(find(uLEDPulses.list_of_durations == duration_pulse))';
    CONDITION = CONDITION_LIST(unique(uLEDPulses.conditionID) == CONDITION_LIST);
      
    uLEDPulses.conditionID =uLEDPulses.conditionID(find(uLEDPulses.conditionID == CONDITION_LIST));
    uLEDPulses.timestamps = uLEDPulses.timestamps(find(uLEDPulses.conditionID == CONDITION_LIST),:);
    uLEDPulses.code = uLEDPulses.code(find(uLEDPulses.conditionID == CONDITION_LIST));
    
    uLEDPulses.list_of_epochs = 1;
    uLEDPulses.list_of_conditions = 1;
    uLEDPulses.list_of_durations(find(uLEDPulses.list_of_durations ~= duration_pulse)) =[];
    % uLEDPulses.conditionID = ones(length(uLEDPulses.conditionID),1); 
    uLEDPulses.conditionID = ones(length(uLEDPulses.conditionID(find(uLEDPulses.conditionID == CONDITION_LIST))),1); 

   % % INDEX = find(round(diff(uLEDPulses.timestamps'),3)==0.02)';
   % uLEDPulses.timestamps = uLEDPulses.timestamps(find(round(diff(uLEDPulses.timestamps'),3)==0.02)',:);
   % uLEDPulses.conditionID =uLEDPulses.conditionID(find(round(diff(uLEDPulses.timestamps'),3)==0.02));
   % uLEDPulses.code = uLEDPulses.code(find(round(diff(uLEDPulses.timestamps'),3)==0.02));
      
end

ints = [];
session = loadSession;
if isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = 'averageCCG_post';
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

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if any(restrict_ints ~= [0 Inf])
    warning('Restricting analysis for intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},restrict_ints);
        spikes.times{ii} = spikes.times{ii}(status);
    end 

    status = InIntervals(uLEDPulses.timestamps(:,1),restrict_ints);
    uLEDPulses.timestamps(~status,:) = [];
    uLEDPulses.conditionID(~status) = [];
    uLEDPulses.code(~status) = [];
end

if length(uLEDPulses.list_of_durations)>length(onset) && length(onset)==1
    onset(1:length(uLEDPulses.list_of_durations)) = onset(1);
end
if length(uLEDPulses.list_of_durations)>length(offset) && length(offset)==1
    offset(1:length(uLEDPulses.list_of_durations)) = offset(1);
end

% remove unncesary fields
spikes_temp = spikes;
clear spikes
spikes.times = spikes_temp.times;
spikes.shankID = spikes_temp.shankID;
spikes.UID = spikes_temp.UID;
clear spikes_temp

codes = 1:max(uLEDPulses.code);

timestamps_recording = min(uLEDPulses.timestamps(:,2)):1/1250:max(uLEDPulses.timestamps(:,2));
if verbose
    disp('Computing cell responses...');          
end

for kk = 1:length(uLEDPulses.list_of_conditions)

    pulseDuration = uLEDPulses.list_of_durations(kk);
    epoch = uLEDPulses.list_of_epochs(kk);
    condition = uLEDPulses.list_of_conditions(kk);

    % generate random events for boostraping
    nPulses = int32(length(find(uLEDPulses.conditionID == uLEDPulses.list_of_conditions(kk)))/...
            length(unique(uLEDPulses.code(uLEDPulses.conditionID == uLEDPulses.list_of_conditions(kk)))));
    randomEvents = [];
    if verbose
        disp('Generating boostrap template...');
    end
    for mm = 1:numRep
        randomEvents{mm} = sort(randsample(timestamps_recording, nPulses))';
    end
    if verbose
        disp('Computing CCG...');
    end
    fprintf('\n');

    [stccg, t] = CCG([spikes.times randomEvents],[],'binSize',binSize,'duration',winSize,'norm','rate');
    fprintf('\n');
    if verbose
        disp('Done!');
    end
   
    t_duringPulse = t > 0 + onset(kk) & t < pulseDuration + offset(kk); 
    randomRatesDuringPulse = squeeze(mean(stccg(t_duringPulse, length(spikes.UID)+1:end,1:length(spikes.UID)),1));
    uLEDResponses_interval.bootsTrapRate(:,kk) = mean(randomRatesDuringPulse,1);
    uLEDResponses_interval.bootsTrapRateStd(:,kk) = std(randomRatesDuringPulse,[],1);
    uLEDResponses_interval.bootsTrapRateSEM(:,kk) = std(randomRatesDuringPulse,[],1)/sqrt(numRep);
    if ~isempty(randomRatesDuringPulse)
        for jj = 1:size(randomRatesDuringPulse,2)
            pd = fitdist(randomRatesDuringPulse(:,jj),'normal');
            uLEDResponses_interval.bootsTrapCI(jj,kk,1:2) = pd.icdf(bootsTrapCI);
        end
    else
        uLEDResponses_interval.bootsTrapCI(1:size(randomRatesDuringPulse,2),kk,1:2) = NaN;
    end
    if verbose
        disp('Collecting responses...');
    end

    for jj = 1:length(codes)
        
        pulses = uLEDPulses.timestamps(uLEDPulses.code == codes(jj) & uLEDPulses.conditionID==kk,1);
        status = InIntervals(pulses, intervals);
        
        % for boostraping, create fake status

        if strcmpi(boostraping_type, 'interval')
        clear fake_status
           nIntervals = size(intervals,1);
           dur_interval = median(diff(intervals'));
           if ~isempty(pulses) 
               fake_train_of_intervals = intervals(1,1):0.001:intervals(end,1);
                fake_status = cell(1, numRep_fake); % Preallocate outside parfor
        
              parfor mm = 1:numRep_fake
                    rand_intervals = sort(fake_train_of_intervals(randperm(length(fake_train_of_intervals), nIntervals)));
                  fake_status{mm} = InIntervals(pulses, [rand_intervals' rand_intervals'+ dur_interval]);            
                 end
             else
                for mm = 1:numRep_fake
                     fake_status{mm} = [0];            
                 end
             end
         elseif strcmpi(boostraping_type, 'pulses')
             clear fake_status
             nPulses = size(pulses,1);
            nStatus = length(find(status==1));
         
             for mm = 1:numRep_fake
                 % disp(mm);
                 rand_status = zeros(nPulses,1);
                 rand_status(randperm(nPulses, nStatus)) = 1;
                 fake_status{mm} = rand_status;            
             end
         
         else
             error('Boostraping type do not recognized! ');
         end
        
          covert fake_status to times
         for ii = 1:length(fake_status)
             fake_status_temp{ii} = pulses(find(fake_status{ii}));
         end
         fake_status = fake_status_temp;
         if isempty(pulses) || length(find(status)) < minNumberOfPulses % 0?
             for iii = 1:numRep_fake
                 fake_status{iii} = [0];
             end
         end
        
        times = spikes.times;        
        if ~isempty(pulses) && length(pulses) > minNumberOfPulses
            times{length(times)+1} = pulses(status==1,1); 
            times{length(times)+1} = pulses(status==0,1); 
            
        else
            times{length(times)+1} = [0]; 
            times{length(times)+1} = [0]; 
        end
        
        % times = cat(2,times, fake_status);
        times = cat(2,times);



        fprintf('\n');
        [stccg, t] = CCG(times,[],'binSize',binSize,'duration',winSize,'norm','rate');
        fprintf('\n');
        numberOfcells = length(spikes.shankID);

        in_interval.responsecurve(:,kk,jj,:) = squeeze(stccg(:,numberOfcells + 1 , 1:numberOfcells))';
        out_interval.responsecurve(:,kk,jj,:) = squeeze(stccg(:, numberOfcells + 2 , 1:numberOfcells))';
        
        if length(times{numberOfcells + 1}) < minNumberOfPulses
            in_interval.responsecurve(:,kk,jj,:) = in_interval.responsecurve(:,kk,jj,:) * NaN;
        end
        if length(times{numberOfcells + 2}) < minNumberOfPulses
            out_interval.responsecurve(:,kk,jj,:) = out_interval.responsecurve(:,kk,jj,:) * NaN;
        end
    %% interpolation in e out
     if interpolate_pulse_sides   
        samples_to_interpolate = [find(t==0)-1:find(t==0)+1; ...
        find(t==pulseDuration)-1:find(t==pulseDuration)+1];

        for ii = 1:size(samples_to_interpolate,1)
            x_axis = 1:size(t,1);
            x_axis(samples_to_interpolate(ii,:)) = [];
            for mm = 1:size(out_interval.responsecurve,1)
                 out_interval.responsecurve(mm,kk,jj,samples_to_interpolate(ii,:)) = ...
                    interp1(x_axis,squeeze(out_interval.responsecurve(mm,1,jj,x_axis)),samples_to_interpolate(ii,:));
               in_interval.responsecurve(mm,1,jj,samples_to_interpolate(ii,:)) = ...
                    interp1(x_axis,squeeze(in_interval.responsecurve(mm,1,jj,x_axis)),samples_to_interpolate(ii,:));
            end
        end
      end
        t_duringPulse = t > 0 + onset(kk) & t < pulseDuration + offset(kk); 
        t_beforePulse = t > -pulseDuration & t < 0; 
        
        in_interval = aux_computeResponse(in_interval, kk, jj, t, t_duringPulse, t_beforePulse, codes, uLEDResponses_interval, pulseDuration, epoch, condition, pulses(status==1,1), spikes, getRaster);
        out_interval = aux_computeResponse(out_interval, kk, jj, t, t_duringPulse, t_beforePulse, codes, uLEDResponses_interval, pulseDuration, epoch, condition, pulses(status==0,1), spikes, getRaster);


        for ii = 1:size(in_interval.responsecurve,1)
            try 
                [h,p]= kstest2(squeeze(in_interval.responsecurve(ii,kk,jj,t_beforePulse)),...
                    squeeze(out_interval.responsecurve(ii,kk,jj,t_beforePulse)));
                h = ~h;
            catch
                h = NaN;
                p = NaN;
            end
            in_interval.is_rateBeforePulse_similar_h(ii,kk,jj,1) = double(h);
            in_interval.is_rateBeforePulse_similar_p(ii,kk,jj,1) = p;
    
            out_interval.is_rateBeforePulse_similar_h(ii,kk,jj,1) = double(h);
            out_interval.is_rateBeforePulse_similar_p(ii,kk,jj,1) = p;
        end
        
        % computation for boostraping
        temp = permute(squeeze(stccg(:, numberOfcells+3:end, 1:numberOfcells)), [3,1,2]);
        
        for ii = 1:numRep
            rand_interval{ii}.responsecurve(:,kk,jj,:) = temp(:,:,ii);
        end
        
        if length(times{numberOfcells + 1}) < minNumberOfPulses
            for ii = 1:numRep
                rand_interval{ii}.responsecurve(:,kk,jj,:) = temp(:,:,ii)*NaN;
            end
        end 

        
       % interpolation rand 
        if interpolate_pulse_sides
            for ii = 1:size(samples_to_interpolate,1)
                x_axis = 1:size(t,1);
                x_axis(samples_to_interpolate(ii,:)) = [];
                for ss=1: numRep
                    for mm = 1:size(out_interval.responsecurve,1)
                         rand_interval{ss}.responsecurve(mm,kk,jj,samples_to_interpolate(ii,:)) = ...
                            interp1(x_axis,squeeze(rand_interval{ss}.responsecurve(mm,kk,jj,x_axis)),samples_to_interpolate(ii,:));
                    end
                end
             end
         end
         for ii = 1:numRep
            rand_interval{ii} = aux_computeResponse(rand_interval{ii}, kk, jj, t, t_duringPulse, t_beforePulse, codes, uLEDResponses_interval, pulseDuration, epoch, condition, times{numberOfcells+2+ii}, spikes, getRaster);
         end 
 
          for ii = 1:numRep
                rand_interval{ii}.is_rateBeforePulse_similar_h(mm,kk,jj,1) = double(h);
                rand_interval{ii}.is_rateBeforePulse_similar_p(mm,kk,jj,1) = p;
           end
        end
end


 %%%%%%%%% SEAL OF MARTA APPORVAL -written by Mario
out_interval.timestamps = t;
in_interval.timestamps = t;
uLEDResponses_interval.timestamps = t;

disp('Parsing cells responses...');

% parse cell responses
for kk = 1:length(uLEDPulses.list_of_conditions)

    for ii = 1:length(spikes.UID)
        %% in and out leds
        in_interval.noRespLEDs.LEDs{ii,kk} = find(in_interval.bootsTrapTest(ii,kk,:) == 0);
        in_interval.respLEDs.LEDs{ii,kk} = find(in_interval.bootsTrapTest(ii,kk,:) == 1);
        in_interval.negRespLEDs.LEDs{ii,kk} = find(in_interval.bootsTrapTest(ii,kk,:)== -1);

        ratio = squeeze(in_interval.rateDuringPulse(ii,kk,:)./in_interval.rateBeforePulse(ii,kk,:));
        ratio(isinf(ratio)) = NaN;
        maxRespLED = find(max(ratio)==ratio,1);
        if length(ratio)==1 % if only 1 light, get everything in maxRespLED
            maxRespLED = 1;
        end
        minRespLED = find(min(ratio)==ratio,1);

        if in_interval.bootsTrapTest(ii,kk,maxRespLED) == 1
            in_interval.maxRespLED.LEDs(ii,kk) = maxRespLED;
        else
            in_interval.maxRespLED.LEDs(ii,kk) = NaN;
        end

        if in_interval.bootsTrapTest(ii,kk,minRespLED) == -1
            in_interval.minRespLED.LEDs(ii,kk) = minRespLED;
        else
            in_interval.minRespLED.LEDs(ii,kk) = NaN;
        end

        out_interval.noRespLEDs.LEDs{ii,kk} = find(out_interval.bootsTrapTest(ii,kk,:) == 0);
        out_interval.respLEDs.LEDs{ii,kk} = find(out_interval.bootsTrapTest(ii,kk,:) == 1);
        out_interval.negRespLEDs.LEDs{ii,kk} = find(out_interval.bootsTrapTest(ii,kk,:)== -1);
        
        ratio = squeeze(out_interval.rateDuringPulse(ii,kk,:)./out_interval.rateBeforePulse(ii,kk,:));
        ratio(isinf(ratio)) = NaN;
        maxRespLED = find(max(ratio)==ratio,1);
        if length(ratio)==1 % if only 1 light, get everything in maxRespLED
            maxRespLED = 1;
        end
        minRespLED = find(min(ratio)==ratio,1);
        
        if out_interval.bootsTrapTest(ii,kk,maxRespLED) == 1
            out_interval.maxRespLED.LEDs(ii,kk) = maxRespLED;
        else
            out_interval.maxRespLED.LEDs(ii,kk) = NaN;
        end
        
        if out_interval.bootsTrapTest(ii,kk,minRespLED) == -1
            out_interval.minRespLED.LEDs(ii,kk) = minRespLED;
        else
            out_interval.minRespLED.LEDs(ii,kk) = NaN;
        end

    %
        in_interval = aux_computeResponse_2(in_interval, kk,spikes,t,out_interval,ii,t_duringPulse); 
        out_interval = aux_computeResponse_2(out_interval, kk, spikes,t,out_interval,ii,t_duringPulse);

         for mm = 1:numRep
         rand_interval{mm} = aux_computeResponse_2(rand_interval{mm}, kk,spikes,t,out_interval,ii,t_duringPulse);
        end
    end
    
end

% check here if the code still works


in_interval.ratioBeforeAfter = in_interval.rateDuringPulse./in_interval.rateBeforePulse;
out_interval.ratioBeforeAfter = out_interval.rateDuringPulse./out_interval.rateBeforePulse;

for ii = 1:numRep
 rand_interval{ii}.ratioBeforeAfter = rand_interval{ii}.rateDuringPulse./rand_interval{ii}.rateBeforePulse;
end

% find cells in non-stimulated shanks
if isfield(spikes,'shankID')
    spikes.shankID = ones(size(spikes.UID));
end
try 
    uLEDResponses_interval.unisInNonStimulatedShanks = (spikes.shankID == uLEDPulses.nonStimulatedShank);
catch
    uLEDResponses_interval.unisInNonStimulatedShanks = [];
end

%%

% parse non-responsive and responsive cells
uLEDResponses_interval.drivenCells = [];
uLEDResponses_interval.stronglyDrivenCells = []; % >2SD


for kk = 1:length(uLEDPulses.list_of_conditions)
    for ii = 1:length(spikes.UID)
        % in
        if any(in_interval.bootsTrapTest(ii,kk,:)==1)
            in_interval.drivenCells(ii,kk) = true;
        else
            in_interval.drivenCells(ii,kk) = false;
        end

        if any(in_interval.bootsTrapTest(ii,kk,:)==-1)
            in_interval.inhibitedCells(ii,kk) = true;
        else
            in_interval.inhibitedCells(ii,kk) = false;
        end
        % out
        if any(out_interval.bootsTrapTest(ii,kk,:)==1)
            out_interval.drivenCells(ii,kk) = true;
        else
            out_interval.drivenCells(ii,kk) = false;
        end

        if any(out_interval.bootsTrapTest(ii,kk,:)==-1)
            out_interval.inhibitedCells(ii,kk) = true;
        else
            out_interval.inhibitedCells(ii,kk) = false;
        end

        % uLEDResponses_interval
        if any(out_interval.bootsTrapTest(ii,kk,:)==1) || any(in_interval.bootsTrapTest(ii,kk,:)==1)
            uLEDResponses_interval.drivenCells(ii,kk) = true;
        else
            uLEDResponses_interval.drivenCells(ii,kk) = false;
        end

        if any(out_interval.bootsTrapTest(ii,kk,:)==-1) || any(in_interval.bootsTrapTest(ii,kk,:)==-1)
            uLEDResponses_interval.inhibitedCells(ii,kk) = true;
        else
            uLEDResponses_interval.inhibitedCells(ii,kk) = false;
        end
    end
end
uLEDResponses_interval.stronglyDrivenCells = out_interval.maxRespLED.rateZ > 2 | in_interval.maxRespLED.rateZ > 2;

in_interval_raster = in_interval;
out_interval_raster = out_interval;
in_interval = rmfield(in_interval,'raster');
out_interval = rmfield(out_interval,'raster');

uLEDResponses_interval_raster = uLEDResponses_interval;
uLEDResponses_interval_raster.in_interval = in_interval_raster;
uLEDResponses_interval_raster.out_interval = out_interval_raster;

uLEDResponses_interval.in_interval = in_interval;
uLEDResponses_interval.out_interval = out_interval;

Temp_rand ={};
for mm= 1 : numRep
    Temp_rand{mm}.rate = rand_interval{mm}.maxRespLED.rate;
    Temp_rand{mm}.rateZ = rand_interval{mm}.maxRespLED.rateZ;
end 

uLEDResponses_interval.rand_interval = Temp_rand;

uLEDResponses_interval.is_rateBeforePulse_similar_h = in_interval.maxRespLED.is_rateBeforePulse_similar_h;
uLEDResponses_interval.is_rateBeforePulse_similar_p = in_interval.maxRespLED.is_rateBeforePulse_similar_p;



% flattening in_interval and out_interval maxLED responses
uLEDResponses_interval.in_interval_maxRespLED = uLEDResponses_interval.in_interval.maxRespLED;
uLEDResponses_interval.out_interval_maxRespLED = uLEDResponses_interval.out_interval.maxRespLED;

if saveMat
    disp('Saving...');
    save([basenameFromBasepath(pwd) '.' save_as '.cellinfo.mat'],'uLEDResponses_interval');
    if getRaster
        save([basenameFromBasepath(pwd) '.' save_as '_raster.cellinfo.mat'],'uLEDResponses_interval_raster','-v7.3');
    end
end

if doPlot
    statisticDots = linspace(-0.015, -0.005, 4);
    nLEDS = size(uLEDResponses_interval.in_interval.responsecurve,3);
    for kk = 1:length(uLEDPulses.list_of_conditions)
        figure;
        set(gcf,'Position',[100 -600 2500 1200]);
        tiledlayout(10,ceil(size(spikes.UID,2)/10),'TileSpacing','tight','Padding','tight');
        for ii = 1:size(uLEDResponses_interval.bootsTrapCI,1)
            % subplot(10,ceil(size(spikes.UID,2)/10),ii); %
            nexttile
            imagesc(uLEDResponses_interval.in_interval.timestamps, 1:nLEDS, squeeze(uLEDResponses_interval.in_interval.responsecurve(ii,kk,:,:))); caxis([0 30]); 
            xlim(winSizePlot); ylim([-0.5 nLEDS+0.5])
            hold on
            plot([0 0], [0 nLEDS+2],'w','LineWidth',1.5);
            plot([uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1) uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1)], [0 nLEDS+2],'color',[.5 .5 .5],'LineWidth',1.5);
            title(num2str(ii),'FontWeight','normal','FontSize',10);
            set(gca,'TickDir','out')
            
            if ii == 1
                ylabel('LED [#][s][0-30Hz]');
            end
            
            modulationSignificanceLevel = squeeze(uLEDResponses_interval.in_interval.modulationSignificanceLevel(ii,kk,:));
            bootsTrapTest = squeeze(uLEDResponses_interval.in_interval.bootsTrapTest(ii,kk,:));
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
        mkdir('SummaryFigures');
        % exportgraphics(gcf,['SummaryFigures\uLEDResponse_rate_condition',num2str(kk),'_dur',num2str(uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1)),'s.png']);
        
        figure;
        set(gcf,'Position',[100 -600 2500 1200]);
        tiledlayout(10,ceil(size(spikes.UID,2)/10),'TileSpacing','tight','Padding','tight');
        for ii = 1:size(uLEDResponses_interval.bootsTrapCI,1)
            % subplot(10,ceil(size(spikes.UID,2)/10),ii); %
            nexttile
            imagesc(uLEDResponses_interval.in_interval.timestamps, 1:nLEDS, squeeze(uLEDResponses_interval.in_interval.responsecurveZ(ii,kk,:,:))); caxis([-5 5]); 
            xlim(winSizePlot); ylim([-0.5 nLEDS+0.5])
            hold on
            plot([0 0], [0 nLEDS+2],'w','LineWidth',1.5);
            plot([uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1) uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1)], [0 nLEDS+2],'color',[.5 .5 .5],'LineWidth',1.5);
            title(num2str(ii),'FontWeight','normal','FontSize',10);
            set(gca,'TickDir','out')
            
            if ii == 1
                ylabel('LED [#][s][+-5 s.d.]');
            end
            
            modulationSignificanceLevel = squeeze(uLEDResponses_interval.in_interval.modulationSignificanceLevel(ii,kk,:));
            bootsTrapTest = squeeze(uLEDResponses_interval.in_interval.bootsTrapTest(ii,kk,:));
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
        % exportgraphics(gcf,['SummaryFigures\uLEDResponse_Zscored_condition',num2str(kk),'_dur',num2str(uLEDResponses_interval.in_interval.pulseDuration(ii,kk,1)),'s.png']);
    end

    % comparison in vs out
    figure
    nexttile
    groupStats({uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1),...
        uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1)},[],'plotType','roundPlot','inAxis',true,'Color',[.2 .2 .8; .2 .6 .6]);
    posData = randn(length(uLEDResponses_interval.in_interval.rateDuringPulse(uLEDResponses_interval.drivenCells==1)),1)/10; 
    posData((posData)>0.3) = posData((posData)>0.3)/2;
    posData((posData)<-0.3) = posData((posData)<-0.3)/2;
    in_interval_rate = uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1);
    out_interval_rate = uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1);

    hold on
    plot(posData+1,uLEDResponses_interval.in_interval.rateDuringPulse(uLEDResponses_interval.drivenCells==1),'o','color',[1 1 1],...
                       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3);
    plot(posData+2,uLEDResponses_interval.out_interval.rateDuringPulse(uLEDResponses_interval.drivenCells==1),'o','color',[1 1 1],...
                       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3);
    for ii = 1:length(in_interval_rate)
        plot([posData(ii)+1 posData(ii)+2], [in_interval_rate(ii) out_interval_rate(ii)],'color',[.5 .5 .5 .5])
    end
    ylim([min([in_interval_rate; out_interval_rate]) max([in_interval_rate; out_interval_rate])]);
    ylabel('Rate (Hz)'); 
    set(gca,'XTick',[1 2],'XTickLabel',{'In interval', 'Out interval'},'XTickLabelRotation',45)
    
    nexttile
    groupStats({uLEDResponses_interval.in_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1),...
        uLEDResponses_interval.out_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1)},[],'plotType','roundPlot','inAxis',true,'Color',[.2 .2 .8; .2 .6 .6]);
    posData = randn(length(uLEDResponses_interval.in_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1)),1)/10; 
    posData((posData)>0.3) = posData((posData)>0.3)/2;
    posData((posData)<-0.3) = posData((posData)<-0.3)/2;
    in_interval_rate = uLEDResponses_interval.in_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1);
    out_interval_rate = uLEDResponses_interval.out_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1);

    hold on
    plot(posData+1,uLEDResponses_interval.in_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1),'o','color',[1 1 1],...
                       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3);
    plot(posData+2,uLEDResponses_interval.out_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1),'o','color',[1 1 1],...
                       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3);
    for ii = 1:length(in_interval_rate)
        plot([posData(ii)+1 posData(ii)+2], [in_interval_rate(ii) out_interval_rate(ii)],'color',[.5 .5 .5 .5])
    end
    ylim([min([in_interval_rate; out_interval_rate]) max([in_interval_rate; out_interval_rate])]);
    ylabel('Rate (Z)'); 
    set(gca,'XTick',[1 2],'XTickLabel',{'In interval', 'Out interval'},'XTickLabelRotation',45);
    
    nexttile
    hold on
    histogram(uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1));
    histogram(uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1));
    legend('Out','In');
    xlabel('Rate (Hz)'); ylabel('#');

    nexttile
    hold on
    histogram(uLEDResponses_interval.out_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1));
    histogram(uLEDResponses_interval.in_interval.maxRespLED.rateZ(uLEDResponses_interval.drivenCells==1));
    legend('Out','In');
    xlabel('Rate (Z)'); ylabel('#');

    nexttile
    hold on
    histogram(uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1) - uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1));
    xlabel('Out - In (Hz)');  ylabel('#');

    nexttile
    hold on
    histogram((uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1) - uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1))./ ...
        (uLEDResponses_interval.out_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1) + uLEDResponses_interval.in_interval.maxRespLED.rate(uLEDResponses_interval.drivenCells==1)));
    xlabel('(Out - In)/ (Out + In) (Hz)');  ylabel('#');

    exportgraphics(gcf,['SummaryFigures\' save_as '_comparison',num2str(kk),'_dur',num2str(uLEDResponses_interval.in_interval.pulseDuration(1,kk,1)),'s.png']);

end
cd(prevPath);

end

%% FIRST FUNCTION

function xx_interval = aux_computeResponse(x_interval, exp_condition, uled_code, t, t_duringPulse, t_beforePulse, codes, uLEDResponses_interval, pulseDuration, epoch, condition, pulses, spikes, getRaster)
    
kk = exp_condition;
jj = uled_code;

for ii = 1:size(x_interval.responsecurve,1)
    % x_interval puo essere in_interval o out_interval
    x_interval.responsecurveZ(ii,kk,jj,:) = (x_interval.responsecurve(ii,kk,jj,:) - mean(x_interval.responsecurve(ii,kk,jj,t < 0)))...
        /std(x_interval.responsecurve(ii,kk,jj,t < 0));
    x_interval.responsecurveZSmooth(ii,kk,jj,:) = smooth(x_interval.responsecurveZ(ii,kk,jj,:));
    x_interval.rateDuringPulse(ii,kk,jj,1) = mean(x_interval.responsecurve(ii,kk,jj,t_duringPulse));
    x_interval.rateBeforePulse(ii,kk,jj,1) = mean(x_interval.responsecurve(ii,kk,jj,t_beforePulse));
    x_interval.rateZDuringPulse(ii,kk,jj,1) = mean(x_interval.responsecurveZ(ii,kk,jj,t_duringPulse));
    x_interval.rateZBeforePulse(ii,kk,jj,1) = mean(x_interval.responsecurveZ(ii,kk,jj,t_beforePulse));
    x_interval.codes(ii,kk,jj,1) = codes(jj);
    
    try 
        [h, x_interval.modulationSignificanceLevel(ii,kk,jj,1)] = kstest2(squeeze(x_interval.responsecurve(ii,kk,jj,t_duringPulse))...
            ,squeeze(x_interval.responsecurve(ii,kk,jj,t_beforePulse)));
    catch
        x_interval.modulationSignificanceLevel(ii,kk,jj,1) = NaN;
    end
    
    ci = squeeze(uLEDResponses_interval.bootsTrapCI(ii,kk,:));
    if x_interval.rateDuringPulse(ii,kk,jj,1) > ci(2)
        test = 1;
    elseif x_interval.rateDuringPulse(ii,kk,jj,1) < ci(1)
        test = -1;
    elseif isnan(x_interval.rateDuringPulse(ii,kk,jj,1))
        test = NaN;
    else
        test = 0;
    end
    x_interval.bootsTrapTest(ii,kk,jj,1) = test;

    if x_interval.rateZDuringPulse(ii,kk,jj,1)  > 1.96
            test = 1;
    elseif x_interval.rateZDuringPulse(ii,kk,jj,1)  < -1.96
        test = -1;
    elseif isnan(x_interval.rateZDuringPulse(ii,kk,jj,1))
        test = NaN;
    else
        test = 0;
    end
    x_interval.zscoreTest(ii,kk,jj,1) = test;
    x_interval.pulseDuration(ii,kk,jj,1) = pulseDuration;
    x_interval.epoch(ii,kk,jj,1) = epoch;
    x_interval.condition(ii,kk,jj,1) = condition;
    x_interval.numberOfPulses(ii,kk,jj,1) = length(pulses);
    
    rasterX = [];
    rasterY = [];
    if getRaster
        for zz = 1:size(pulses,1)
            temp_spk = spikes.times{ii}(find(spikes.times{ii} - pulses(zz,1)  > salt_time(1) & spikes.times{ii} - pulses(zz,1)  < salt_time(2))) - pulses(zz,1);
            rasterX = [rasterX; temp_spk];
            if ~isempty(temp_spk)
                rasterY = [rasterY; zz * ones(size((temp_spk)))];
            end
        end
    end

    if ~isempty(rasterX)
        [rasterHist3,c] = hist3([rasterY rasterX],{1:size(pulses,1) salt_time(1):salt_binSize:salt_time(2)});
        x_interval.raster.rasterCount{ii,kk,jj} = rasterHist3;
        x_interval.raster.rasterProb{ii,kk,jj} = rasterHist3/sum(rasterHist3(:));
        x_interval.raster.TrialsNumber{ii,kk,jj} = c{1};
        x_interval.raster.times{ii,kk,jj} = c{2};
        x_interval.raster.rasterTrials{ii,kk,jj} = rasterY;
        x_interval.raster.rasterSpikesTimes{ii,kk,jj} = rasterX;

        time = c{2};
        baseidx = dsearchn(time', salt_baseline');
        tidx = dsearchn(time', [0; pulseDuration*2]);    
        st = length(baseidx(1):baseidx(2));
        nmbn = round(salt_win/salt_binSize);
        v = 1:nmbn:st;
        if any((v + nmbn - 1) > st)
            error('reduce window size or baseline duration')
        end
        [x_interval.salt.p_value(ii,kk,jj,1), x_interval.salt.I_statistics(ii,kk,jj,1)] = salt(rasterHist3(:,baseidx(1):baseidx(2)),rasterHist3(:,tidx(1):tidx(2)),salt_binSize, salt_win);
    else
        x_interval.raster.rasterCount{ii,kk,jj} = NaN;
        x_interval.raster.rasterProb{ii,kk,jj} = NaN;
        x_interval.raster.TrialsNumber{ii,kk,jj} = NaN;
        x_interval.raster.times{ii,kk,jj} = NaN;
        x_interval.raster.rasterTrials{ii,kk,jj} = NaN;
        x_interval.raster.rasterSpikesTimes{ii,kk,jj} = NaN;
        x_interval.salt.p_value(ii,kk,jj,1) = NaN;
        x_interval.salt.I_statistics(ii,kk,jj,1) = NaN;
    end
end
xx_interval = x_interval;

end
%% SECOND FUNCTION 

function x_interval = aux_computeResponse_2(x_interval,kk,spikes,t,interval_out,ii,t_duringPulse)

out_interval=interval_out;

% for ii = 1:length(spikes.UID)
    % noRespLEDs
    x_interval.noRespLEDs.rate{ii,kk} = x_interval.rateDuringPulse(ii,kk,out_interval.noRespLEDs.LEDs{ii,kk});
    x_interval.noRespLEDs.rateBeforePulse{ii,kk} = x_interval.rateBeforePulse(ii,kk,out_interval.noRespLEDs.LEDs{ii,kk});
    x_interval.noRespLEDs.rateZ{ii,kk} = x_interval.rateZDuringPulse(ii,kk,out_interval.noRespLEDs.LEDs{ii,kk});
    x_interval.noRespLEDs.rateZBeforePulse{ii,kk} = x_interval.rateZBeforePulse(ii,kk,out_interval.noRespLEDs.LEDs{ii,kk});
    x_interval.noRespLEDs.ratioBeforeAfter{ii,kk} = x_interval.noRespLEDs.rate{ii,kk}./x_interval.noRespLEDs.rateBeforePulse{ii,kk}; 
    x_interval.noRespLEDs.meanRatio(ii,kk) = nanmean(x_interval.noRespLEDs.ratioBeforeAfter{ii,kk});
    x_interval.noRespLEDs.meanRateBeforePulse(ii,kk) = nanmean(x_interval.noRespLEDs.rateBeforePulse{ii,kk});
    x_interval.noRespLEDs.meanRateDuringPulse(ii,kk) = nanmean(x_interval.noRespLEDs.rate{ii,kk});
    if isnan(x_interval.noRespLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
        x_interval.noRespLEDs.LEDs{ii,kk} = [];
        x_interval.noRespLEDs.rate{ii,kk} = [];
        x_interval.noRespLEDs.rateBeforePulse{ii,kk} = [];
        x_interval.noRespLEDs.rateZ{ii,kk} = [];
        x_interval.noRespLEDs.rateZBeforePulse{ii,kk} = [];
        x_interval.noRespLEDs.ratioBeforeAfter{ii,kk} = [];
    end

    % negRespLEDs
    x_interval.negRespLEDs.rate{ii,kk} = x_interval.rateDuringPulse(ii,kk,out_interval.negRespLEDs.LEDs{ii,kk});
    x_interval.negRespLEDs.rateBeforePulse{ii,kk} = x_interval.rateBeforePulse(ii,kk,out_interval.negRespLEDs.LEDs{ii,kk});  
    x_interval.negRespLEDs.rateZ{ii,kk} = x_interval.rateZDuringPulse(ii,kk,out_interval.negRespLEDs.LEDs{ii,kk});
    x_interval.negRespLEDs.rateZBeforePulse{ii,kk} = x_interval.rateZBeforePulse(ii,kk,out_interval.negRespLEDs.LEDs{ii,kk});  
    x_interval.negRespLEDs.ratioBeforeAfter{ii,kk} = x_interval.negRespLEDs.rate{ii,kk}./x_interval.negRespLEDs.rateBeforePulse{ii,kk};  
    x_interval.negRespLEDs.meanRatio(ii,kk) = nanmean(x_interval.negRespLEDs.ratioBeforeAfter{ii,kk});
    x_interval.negRespLEDs.ratioNoResp(ii,kk) = nanmean(x_interval.negRespLEDs.rate{ii,kk})./nanmean(x_interval.noRespLEDs.rate{ii,kk});
    x_interval.negRespLEDs.meanRateBeforePulse(ii,kk) = nanmean(x_interval.negRespLEDs.rateBeforePulse{ii,kk});
    x_interval.negRespLEDs.meanRateDuringPulse(ii,kk) = nanmean(x_interval.negRespLEDs.rate{ii,kk});
    if isnan(x_interval.negRespLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
        x_interval.negRespLEDs.LEDs{ii,kk} = [];
        x_interval.negRespLEDs.rate{ii,kk} = [];
        x_interval.negRespLEDs.rateBeforePulse{ii,kk} = [];
        x_interval.negRespLEDs.rateZ{ii,kk} = [];
        x_interval.negRespLEDs.rateZBeforePulse{ii,kk} = [];
        x_interval.negRespLEDs.ratioBeforeAfter{ii,kk} = [];
    end

    % respLEDs
    x_interval.respLEDs.rate{ii,kk} = x_interval.rateDuringPulse(ii,kk,out_interval.respLEDs.LEDs{ii,kk});
    x_interval.respLEDs.rateBeforePulse{ii,kk} = x_interval.rateBeforePulse(ii,kk,out_interval.respLEDs.LEDs{ii,kk});  
    x_interval.respLEDs.rateZ{ii,kk} = x_interval.rateZDuringPulse(ii,kk,out_interval.respLEDs.LEDs{ii,kk});
    x_interval.respLEDs.rateZBeforePulse{ii,kk} = x_interval.rateZBeforePulse(ii,kk,out_interval.respLEDs.LEDs{ii,kk});  
    x_interval.respLEDs.ratioBeforeAfter{ii,kk} = x_interval.respLEDs.rate{ii,kk}./x_interval.respLEDs.rateBeforePulse{ii,kk};  
    x_interval.respLEDs.meanRatio(ii,kk) = nanmean(x_interval.respLEDs.ratioBeforeAfter{ii,kk});
    x_interval.respLEDs.ratioNoResp(ii,kk) = nanmean(x_interval.respLEDs.rate{ii,kk})./nanmean(x_interval.noRespLEDs.rate{ii,kk});
    x_interval.respLEDs.meanRateBeforePulse(ii,kk) = nanmean(x_interval.respLEDs.rateBeforePulse{ii,kk});
    x_interval.respLEDs.meanRateDuringPulse(ii,kk) = nanmean(x_interval.respLEDs.rate{ii,kk});
   
    %%
    if length(out_interval.respLEDs.LEDs{ii,kk}) == 1
        dim_mean = 2;
    else
        dim_mean = 1;
    end
    x_interval.respLEDs.responseCurve(ii,kk,:) = nanmean(squeeze(x_interval.responsecurve(ii,kk,out_interval.respLEDs.LEDs{ii,kk},:)),dim_mean);
    x_interval.respLEDs.responseCurveZ(ii,kk,:) = nanmean(squeeze(x_interval.responsecurveZ(ii,kk,out_interval.respLEDs.LEDs{ii,kk},:)),dim_mean);
    x_interval.respLEDs.responseCurveZSmooth(ii,kk,:) = nanmean(squeeze(x_interval.responsecurveZSmooth(ii,kk,out_interval.respLEDs.LEDs{ii,kk},:)),dim_mean);
    if isnan(x_interval.respLEDs.meanRatio(ii,kk)) % if only NaN, remove LEDs
        x_interval.respLEDs.LEDs{ii,kk} = [];
        x_interval.respLEDs.rate{ii,kk} = [];
        x_interval.respLEDs.rateBeforePulse{ii,kk} = [];
        x_interval.respLEDs.rateZ{ii,kk} = [];
        x_interval.respLEDs.rateZBeforePulse{ii,kk} = [];
        x_interval.respLEDs.ratioBeforeAfter{ii,kk} = [];
        x_interval.respLEDs.responseCurve(ii,kk,:) = nan(size(t));
        x_interval.respLEDs.responseCurveZ(ii,kk,:) = nan(size(t));
        x_interval.respLEDs.responseCurveZSmooth(ii,kk,:) = nan(size(t));
    end
    
    % maxRespLED % the reference is the out_intervals max leds!!!
    if ~isnan(out_interval.maxRespLED.LEDs(ii,kk))
        x_interval.maxRespLED.rate(ii,kk) = x_interval.rateDuringPulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));
        x_interval.maxRespLED.rateBeforePulse(ii,kk) = x_interval.rateBeforePulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));   
        x_interval.maxRespLED.rateZ(ii,kk) = x_interval.rateZDuringPulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));
        x_interval.maxRespLED.rateZBeforePulse(ii,kk) = x_interval.rateZBeforePulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));   
        x_interval.maxRespLED.ratioBeforeAfter(ii,kk) = x_interval.maxRespLED.rate(ii,kk)./x_interval.maxRespLED.rateBeforePulse(ii,kk); 
        x_interval.maxRespLED.ratioNoResp(ii,kk) = x_interval.is_rateBeforePulse_similar_h(ii,kk);
        
        x_interval.maxRespLED.is_rateBeforePulse_similar_h(ii,kk) = x_interval.is_rateBeforePulse_similar_h(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));
        x_interval.maxRespLED.is_rateBeforePulse_similar_p(ii,kk) = x_interval.is_rateBeforePulse_similar_p(ii,kk,out_interval.maxRespLED.LEDs(ii,kk));

        if isnan(x_interval.maxRespLED.rate(ii,kk))
            x_interval.maxRespLED.LEDs(ii,kk) = NaN;
        end
        x_interval.maxRespLED.responseCurve(ii,kk,:) = squeeze(x_interval.responsecurve(ii,kk,out_interval.maxRespLED.LEDs(ii,kk),:));
        x_interval.maxRespLED.responseCurveZ(ii,kk,:) = squeeze(x_interval.responsecurveZ(ii,kk,out_interval.maxRespLED.LEDs(ii,kk),:));
        x_interval.maxRespLED.responseCurveZSmooth(ii,kk,:) = squeeze(x_interval.responsecurveZSmooth(ii,kk,out_interval.maxRespLED.LEDs(ii,kk),:));
        if mean(x_interval.maxRespLED.responseCurve(ii,kk,t_duringPulse)) ~= x_interval.rateDuringPulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk)) && ~isnan(x_interval.rateDuringPulse(ii,kk,out_interval.maxRespLED.LEDs(ii,kk)))
            warning('Something terrible happened with indexing!!')
            % keyboard;
        end
    else
        x_interval.maxRespLED.rate(ii,kk) = NaN;
        x_interval.maxRespLED.rateBeforePulse(ii,kk) = NaN;
        x_interval.maxRespLED.rateZ(ii,kk) = NaN;
        x_interval.maxRespLED.rateZBeforePulse(ii,kk) = NaN;
        x_interval.maxRespLED.ratioBeforeAfter(ii,kk) = NaN;
        x_interval.maxRespLED.ratioNoResp(ii,kk) = NaN;
        x_interval.maxRespLED.is_rateBeforePulse_similar_h(ii,kk) = NaN;
        x_interval.maxRespLED.is_rateBeforePulse_similar_p(ii,kk) = NaN;
        x_interval.maxRespLED.responseCurve(ii,kk,:) = nan(size(t));
        x_interval.maxRespLED.responseCurveZ(ii,kk,:) = nan(size(t));
        x_interval.maxRespLED.responseCurveZSmooth(ii,kk,:) = nan(size(t));
    end
    
    % minRespLED
    if ~isnan(out_interval.minRespLED.LEDs(ii,kk))
        x_interval.minRespLED.rate(ii,kk) = x_interval.rateDuringPulse(ii,kk,out_interval.minRespLED.LEDs(ii,kk));
        x_interval.minRespLED.rateBeforePulse(ii,kk) = x_interval.rateBeforePulse(ii,kk,out_interval.minRespLED.LEDs(ii,kk));   
        x_interval.minRespLED.rateZ(ii,kk) = x_interval.rateZDuringPulse(ii,kk,out_interval.minRespLED.LEDs(ii,kk));
        x_interval.minRespLED.rateZBeforePulse(ii,kk) = x_interval.rateZBeforePulse(ii,kk,out_interval.minRespLED.LEDs(ii,kk));   
        x_interval.minRespLED.ratioBeforeAfter(ii,kk) = x_interval.minRespLED.rate(ii,kk)./x_interval.minRespLED.rateBeforePulse(ii,kk);  
        x_interval.minRespLED.ratioNoResp(ii,kk) = x_interval.minRespLED.rate(ii,kk)./nanmean(x_interval.noRespLEDs.rate{ii,kk});

        x_interval.minRespLED.is_rateBeforePulse_similar_h(ii,kk) = x_interval.is_rateBeforePulse_similar_h(ii,kk,out_interval.minRespLED.LEDs(ii,kk));
        x_interval.minRespLED.is_rateBeforePulse_similar_p(ii,kk) = x_interval.is_rateBeforePulse_similar_p(ii,kk,out_interval.minRespLED.LEDs(ii,kk));

        if isnan(x_interval.minRespLED.rate(ii,kk))
            x_interval.minRespLED.LEDs(ii,kk) = NaN;
        end
        x_interval.minRespLED.responseCurve(ii,kk,:) = squeeze(x_interval.responsecurve(ii,kk,out_interval.minRespLED.LEDs(ii,kk),:));
        x_interval.minRespLED.responseCurveZ(ii,kk,:) = squeeze(x_interval.responsecurveZ(ii,kk,out_interval.minRespLED.LEDs(ii,kk),:));
        x_interval.minRespLED.responseCurveZSmooth(ii,kk,:) = squeeze(x_interval.responsecurveZSmooth(ii,kk,out_interval.minRespLED.LEDs(ii,kk),:));
    else
        x_interval.minRespLED.rate(ii,kk) = NaN;
        x_interval.minRespLED.rateBeforePulse(ii,kk) = NaN;   
        x_interval.minRespLED.rateZ(ii,kk) = NaN;
        x_interval.minRespLED.rateZBeforePulse(ii,kk) = NaN;   
        x_interval.minRespLED.ratioBeforeAfter(ii,kk) = NaN;  
        x_interval.minRespLED.ratioNoResp(ii,kk) = NaN;
        x_interval.minRespLED.is_rateBeforePulse_similar_h(ii,kk) = NaN;
        x_interval.minRespLED.is_rateBeforePulse_similar_p(ii,kk) = NaN;
        x_interval.minRespLED.responseCurve(ii,kk,:) = nan(size(t));
        x_interval.minRespLED.responseCurveZ(ii,kk,:) = nan(size(t));
        x_interval.minRespLED.responseCurveZSmooth(ii,kk,:) = nan(size(t));
    end
  
end

