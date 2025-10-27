
function [collision_metrics] = get_light_spike_CollisionMetrics(uLEDResponses_interval, varargin)
% [collisionMetrics] = get_light_spike_CollisionMetrics(uLEDResponses_interval, varargin)
%
% Manu, Neural Computation Lab, 2023

% Parse options
p = inputParser;
addRequired(p,'uLEDResponses_interval',@iscell);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'light_response_threshold_factor',2,@isscalar); % only cells that change their rate at least by this factor are considered light responsive
addParameter(p,'light_response_only_pyr',true,@islogical);
addParameter(p,'before_pulse_factor',1.5,@isscalar); % which cells show similar before pulse activity in and out interval? as measured by std around mean
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'label',[]); % string to add to the title for figures and mat file
addParameter(p,'saveMat',true,@islogical); % 
addParameter(p,'rate_change_threshold',7,@isnumeric); % 
addParameter(p,'spikes',[],@isstruct); % 
addParameter(p,'uLEDPulses',getuLEDPulses,@isstruct); % 
addParameter(p,'interpolate_pulse_sides',false,@islogical); % 
addParameter(p,'update_cell_metrics',true,@islogical); % 
addParameter(p,'save_as','lightSpikeCollisions',@ischar);
addParameter(p,'rand_analysis',false,@islogical); % if you put yes you would save also rand interval values and boostatrapping values



parse(p, uLEDResponses_interval, varargin{:});
basepath = p.Results.basepath;
cell_metrics = p.Results.cell_metrics;
uLEDPulses = p.Results.uLEDPulses;
light_response_threshold_factor = p.Results.light_response_threshold_factor;
light_response_only_pyr = p.Results.light_response_only_pyr;
before_pulse_factor = p.Results.before_pulse_factor;
doPlot = p.Results.doPlot;
label = p.Results.label;
saveMat = p.Results.saveMat;
rate_change_threshold = p.Results.rate_change_threshold;
spikes = p.Results.spikes;
interpolate_pulse_sides = p.Results.interpolate_pulse_sides;
update_cell_metrics = p.Results.update_cell_metrics;
save_as = p.Results.save_as;
rand_analysis = p.Results.rand_analysis;


% Deal with inputs
prevPath = pwd;
cd(basepath);

if isempty(cell_metrics)
    cell_metrics = loadCellMetrics;
end

if ~isstring(label)
   label = string(label);
end

if isempty(spikes)
    try
        spikes = loadSpikes;
    catch
        warning('Spikes structure not found! Spiking waveform and CCG...');
    end
end

% if ~isstruct(uLEDPulses) && isnan(uLEDPulses)
%     uLEDPulses = getuLEDPulses;
% end

% stacking data
uLEDResponses_InInterval.presynapticID = [];
uLEDResponses_InInterval.presynapticCellType = [];
uLEDResponses_InInterval.postsynapticID = [];
uLEDResponses_InInterval.maxRateBeforePulse = [];
uLEDResponses_InInterval.maxRatePulse = [];
uLEDResponses_InInterval.maxZBeforePulse = [];
uLEDResponses_InInterval.maxZPulse = [];
uLEDResponses_InInterval.putativeCellType = [];
uLEDResponses_InInterval.is_out_rateBeforePuse_similar = [];
uLEDResponses_InInterval.responsecurve = [];
uLEDResponses_InInterval.responsecurveZ = [];
uLEDResponses_InInterval.responsecurveZSmooth = [];

uLEDResponses_OutInterval.presynapticID = [];
uLEDResponses_OutInterval.presynapticCellType = [];
uLEDResponses_OutInterval.postsynapticID = [];
uLEDResponses_OutInterval.maxRateBeforePulse = [];
uLEDResponses_OutInterval.maxRatePulse = [];
uLEDResponses_OutInterval.maxZBeforePulse = [];
uLEDResponses_OutInterval.maxZPulse = [];
uLEDResponses_OutInterval.putativeCellType = [];
uLEDResponses_OutInterval.is_in_rateBeforePuse_similar = [];
uLEDResponses_OutInterval.responsecurve = [];
uLEDResponses_OutInterval.responsecurveZ = [];
uLEDResponses_OutInterval.responsecurveZSmooth = [];
if rand_analysis
    uLEDResponses_RandInterval.maxRatePulse = [];
    uLEDResponses_RandInterval.maxZPulse = [];
end 

for ii = 1:length(uLEDResponses_interval)
    if strcmpi(cell_metrics.putativeCellType{ii},'Pyramidal Cell')
        cellType = 1;
    elseif strcmpi(cell_metrics.putativeCellType{ii},'Wide Interneuron')
        cellType = 3;
    elseif strcmpi(cell_metrics.putativeCellType{ii},'Narrow Interneuron')
        cellType = 2;
    else
        cellType = 4;
    end
    % in
    uLEDResponses_InInterval.presynapticID = [uLEDResponses_InInterval.presynapticID; ...
        ones(size(uLEDResponses_interval{ii}.bootsTrapRate))*ii];
    uLEDResponses_InInterval.presynapticCellType = [uLEDResponses_InInterval.presynapticCellType; ...
        ones(size(uLEDResponses_interval{ii}.bootsTrapRate))*cellType];
    uLEDResponses_InInterval.maxRateBeforePulse = [uLEDResponses_InInterval.maxRateBeforePulse; ...
        uLEDResponses_interval{ii}.in_interval.maxRespLED.rateBeforePulse];
    uLEDResponses_InInterval.maxRatePulse = [uLEDResponses_InInterval.maxRatePulse; ...
        uLEDResponses_interval{ii}.in_interval.maxRespLED.rate];
    uLEDResponses_InInterval.maxZBeforePulse = [uLEDResponses_InInterval.maxZBeforePulse; ...
        uLEDResponses_interval{ii}.in_interval.maxRespLED.rateZBeforePulse];
    uLEDResponses_InInterval.maxZPulse = [uLEDResponses_InInterval.maxZPulse; ...
        uLEDResponses_interval{ii}.in_interval.maxRespLED.rateZ];
    uLEDResponses_InInterval.putativeCellType = [uLEDResponses_InInterval.putativeCellType; ...
        cell_metrics.putativeCellType'];
    uLEDResponses_InInterval.is_out_rateBeforePuse_similar = [uLEDResponses_InInterval.is_out_rateBeforePuse_similar;...
        uLEDResponses_interval{ii}.is_rateBeforePulse_similar_h];
    uLEDResponses_InInterval.postsynapticID = [uLEDResponses_InInterval.postsynapticID; ...
        (1:size(uLEDResponses_interval{ii}.bootsTrapRate,1))'];
    uLEDResponses_InInterval.responsecurve = [uLEDResponses_InInterval.responsecurve; ...
        squeeze(uLEDResponses_interval{ii}.in_interval.maxRespLED.responseCurve)];

    % out
    uLEDResponses_OutInterval.presynapticID = [uLEDResponses_OutInterval.presynapticID; ...
        ones(size(uLEDResponses_interval{ii}.bootsTrapRate))*ii];
    uLEDResponses_OutInterval.presynapticCellType = [uLEDResponses_OutInterval.presynapticCellType; ...
        ones(size(uLEDResponses_interval{ii}.bootsTrapRate))*cellType];
    uLEDResponses_OutInterval.maxRateBeforePulse = [uLEDResponses_OutInterval.maxRateBeforePulse; ...
        uLEDResponses_interval{ii}.out_interval.maxRespLED.rateBeforePulse];
    uLEDResponses_OutInterval.maxRatePulse = [uLEDResponses_OutInterval.maxRatePulse; ...
        uLEDResponses_interval{ii}.out_interval.maxRespLED.rate];
    uLEDResponses_OutInterval.maxZBeforePulse = [uLEDResponses_OutInterval.maxZBeforePulse; ...
        uLEDResponses_interval{ii}.out_interval.maxRespLED.rateZBeforePulse];
    uLEDResponses_OutInterval.maxZPulse = [uLEDResponses_OutInterval.maxZPulse; ...
        uLEDResponses_interval{ii}.out_interval.maxRespLED.rateZ];
    uLEDResponses_OutInterval.putativeCellType = [uLEDResponses_OutInterval.putativeCellType; ...
        cell_metrics.putativeCellType'];
    uLEDResponses_OutInterval.is_in_rateBeforePuse_similar = [uLEDResponses_OutInterval.is_in_rateBeforePuse_similar;...
        uLEDResponses_interval{ii}.is_rateBeforePulse_similar_h];
    uLEDResponses_OutInterval.postsynapticID = [uLEDResponses_OutInterval.postsynapticID; ...
        (1:size(uLEDResponses_interval{ii}.bootsTrapRate,1))'];
    uLEDResponses_OutInterval.responsecurve = [uLEDResponses_OutInterval.responsecurve; ...
        squeeze(uLEDResponses_interval{ii}.out_interval.maxRespLED.responseCurve)];

 if rand_analysis

    %rand
    temp_rand_rate = [];
    temp_rand_rateZ = [];
    for jj = 1:length(uLEDResponses_interval{ii}.rand_interval)
        temp_rand_rate(:,jj) = uLEDResponses_interval{ii}.rand_interval{jj}.rate;
        temp_rand_rateZ(:,jj) = uLEDResponses_interval{ii}.rand_interval{jj}.rateZ;
    end
    uLEDResponses_RandInterval.maxRatePulse = [uLEDResponses_RandInterval.maxRatePulse; temp_rand_rate];
    uLEDResponses_RandInterval.maxZPulse = [uLEDResponses_RandInterval.maxZPulse; temp_rand_rateZ];
 end
end
 
timestamps = uLEDResponses_interval{1}.in_interval.timestamps;
uLEDResponses_OutInterval.timestamps = timestamps';
uLEDResponses_InInterval.timestamps = timestamps';

t_duringPulse = timestamps > 0 & timestamps <0.02;
uLEDResponses_OutInterval.maxRatePulse = mean(uLEDResponses_OutInterval.responsecurve(:,t_duringPulse),2);
uLEDResponses_InInterval.maxRatePulse = mean(uLEDResponses_InInterval.responsecurve(:,t_duringPulse),2);


if interpolate_pulse_sides
    timestamps = uLEDResponses_interval{1}.in_interval.timestamps;
    pulse_duration = uLEDResponses_interval{1}.in_interval.pulseDuration(1);
    samples_to_interpolate = [find(timestamps==0)-1:find(timestamps==0)+1; ...
        find(timestamps==pulse_duration)-1:find(timestamps==pulse_duration)+1];
    for ii = 1:size(samples_to_interpolate,1)
        x_axis = 1:size(timestamps,1);
        x_axis(samples_to_interpolate(ii,:)) = [];
        for jj = 1:size(uLEDResponses_OutInterval.responsecurve,1)
            uLEDResponses_OutInterval.responsecurve(jj,samples_to_interpolate(ii,:)) = ...
                interp1(x_axis,uLEDResponses_OutInterval.responsecurve(jj,x_axis),samples_to_interpolate(ii,:));
            uLEDResponses_InInterval.responsecurve(jj,samples_to_interpolate(ii,:)) = ...
                interp1(x_axis,uLEDResponses_InInterval.responsecurve(jj,x_axis),samples_to_interpolate(ii,:));
        end
    end
    % the interpolation may change a bit the max responses... recomputing
    timestamps = uLEDResponses_interval{1}.in_interval.timestamps;
    t_duringPulse = timestamps > 0 & timestamps <0.02;
    uLEDResponses_OutInterval.maxRatePulse = mean(uLEDResponses_OutInterval.responsecurve(:,t_duringPulse),2);
    uLEDResponses_InInterval.maxRatePulse = mean(uLEDResponses_InInterval.responsecurve(:,t_duringPulse),2);
end

% z_scoring properly
timestamps = uLEDResponses_interval{1}.in_interval.timestamps;
pulse_duration = uLEDResponses_interval{1}.in_interval.pulseDuration(1);
z_win = InIntervals(timestamps, [min(timestamps) min(timestamps)+ pulse_duration]);
uLEDResponses_InInterval.responsecurveZ = zscore_win(uLEDResponses_InInterval.responsecurve,z_win)';
uLEDResponses_OutInterval.responsecurveZ = zscore_win(uLEDResponses_OutInterval.responsecurve,z_win)';
uLEDResponses_OutInterval.maxZPulse = mean(uLEDResponses_OutInterval.responsecurveZ(:,t_duringPulse),2);
uLEDResponses_InInterval.maxZPulse = mean(uLEDResponses_InInterval.responsecurveZ(:,t_duringPulse),2);

% light responsive cells
uLEDResponses_OutInterval.lightResponsive = ...
    uLEDResponses_OutInterval.maxRatePulse./uLEDResponses_OutInterval.maxRateBeforePulse > light_response_threshold_factor;
uLEDResponses_InInterval.lightResponsive  = ...
    uLEDResponses_InInterval.maxRatePulse./uLEDResponses_InInterval.maxRateBeforePulse > light_response_threshold_factor;
if light_response_only_pyr
    uLEDResponses_OutInterval.lightResponsive = uLEDResponses_OutInterval.lightResponsive & strcmpi(uLEDResponses_OutInterval.putativeCellType,'Pyramidal Cell');
    uLEDResponses_InInterval.lightResponsive  = uLEDResponses_InInterval.lightResponsive  & strcmpi(uLEDResponses_InInterval.putativeCellType,'Pyramidal Cell');
end
% presynaptic cells is interneuron?
uLEDResponses_InInterval.is_pre_interneuron = ...
    uLEDResponses_InInterval.presynapticCellType  == 2 | uLEDResponses_InInterval.presynapticCellType  == 3;
uLEDResponses_OutInterval.is_pre_interneuron = ...
    uLEDResponses_OutInterval.presynapticCellType == 2 | uLEDResponses_OutInterval.presynapticCellType == 3;

%% collision metrics
% metadata
collision_metrics.label = label;
collision_metrics.is_rate_before_similar_ktest = uLEDResponses_OutInterval.is_in_rateBeforePuse_similar;
change_factor = max([uLEDResponses_OutInterval.maxRateBeforePulse,uLEDResponses_InInterval.maxRateBeforePulse]')'...
    ./min([uLEDResponses_OutInterval.maxRateBeforePulse,uLEDResponses_InInterval.maxRateBeforePulse]')';
collision_metrics.change_factor = change_factor;
collision_metrics.is_rate_before_similar_factor = change_factor<before_pulse_factor;
collision_metrics.range = range([uLEDResponses_OutInterval.maxRateBeforePulse,uLEDResponses_InInterval.maxRateBeforePulse]')';
collision_metrics.is_pre_interneuron = uLEDResponses_InInterval.is_pre_interneuron;
collision_metrics.is_lightResponsive = uLEDResponses_OutInterval.lightResponsive;
collision_metrics.putativeCellType = uLEDResponses_InInterval.putativeCellType;
collision_metrics.is_post_pyramidalCell = strcmpi(collision_metrics.putativeCellType,'Pyramidal Cell');

not_same_cell = uLEDResponses_OutInterval.postsynapticID~=uLEDResponses_OutInterval.presynapticCellType;
collision_metrics.candidate_int_pyr_pairs = collision_metrics.is_pre_interneuron==1 & ...
    collision_metrics.is_post_pyramidalCell & not_same_cell;
        % collision_metrics.is_lightResponsive & collision_metrics.is_rate_before_similar_factor & ...
collision_metrics.candidate_pyr_pyr_pairs = collision_metrics.is_pre_interneuron==0 & ...
        collision_metrics.is_post_pyramidalCell & not_same_cell;
        % collision_metrics.is_lightResponsive & collision_metrics.is_rate_before_similar_factor & ...

pre_post_CCG = [];
pre_waveforms = [];
post_waveforms = [];
pre_firingRate = [];
post_firingRate = [];
if ~isempty(spikes)
    [ccg,t] = CCG(spikes.times,[],'binSize',0.0004,'duration',0.12,'norm','rate');
    fprintf('\n');
    for ii = 1:spikes.numcells
        pre_post_CCG    = [pre_post_CCG; squeeze(ccg(:,ii,:))'];
        pre_firingRate  = [pre_firingRate;...
            ones(size(spikes.times'))*length(spikes.times{ii})/(max(spikes.times{ii})-min(spikes.times{ii}))];
        post_firingRate = [post_firingRate; cellfun(@length,spikes.times)./(cellfun(@max,spikes.times)-cellfun(@min,spikes.times))];
        pre_waveforms   = [pre_waveforms; (spikes.filtWaveform{ii}'* ones(size(spikes.times)))'];
        post_waveforms  = [post_waveforms; cat(1,spikes.filtWaveform{:})];
    end
end
collision_metrics.pre_post_CCG    = pre_post_CCG;
collision_metrics.pre_post_CCG_timestamps    = t;
collision_metrics.pre_waveforms   = pre_waveforms;
collision_metrics.post_waveforms  = post_waveforms;
collision_metrics.pre_firingRate  = pre_firingRate;
collision_metrics.post_firingRate = post_firingRate;

% CCG in baseline

pre_post_baseline_CCG = [];
t_limit = uLEDPulses.timestamps(1);
spikes_times_baseline = {};
if ~isempty(spikes)
    for ii = 1:length(spikes.times)
        spikes_times_baseline{ii} = spikes.times{ii}(spikes.times{ii} < t_limit);
    end
    [ccg,t] = CCG(spikes_times_baseline,[],'binSize',0.0004,'duration',0.12,'norm','rate');
    fprintf('\n');
    for ii = 1:spikes.numcells
        pre_post_baseline_CCG    = [pre_post_baseline_CCG; squeeze(ccg(:,ii,:))'];
    end
end
collision_metrics.pre_post_baseline_CCG  = pre_post_baseline_CCG;
% actual metrics
collision_metrics.rate_difference            = uLEDResponses_OutInterval.maxRatePulse - uLEDResponses_InInterval.maxRatePulse;
collision_metrics.rateZ_difference           = uLEDResponses_OutInterval.maxZPulse - uLEDResponses_InInterval.maxZPulse;
collision_metrics.rate_only_light            = uLEDResponses_OutInterval.maxRatePulse;
collision_metrics.rateZ_only_light           = uLEDResponses_OutInterval.maxZPulse;
collision_metrics.rate_before_light          = uLEDResponses_OutInterval.maxRateBeforePulse;
collision_metrics.rateZ_before_light         = uLEDResponses_OutInterval.maxZBeforePulse;
collision_metrics.rate_ligh_spike_collision  = uLEDResponses_InInterval.maxRatePulse;
collision_metrics.rateZ_ligh_spike_collision = uLEDResponses_InInterval.maxZPulse;
% 
preInt_select = collision_metrics.candidate_int_pyr_pairs;
collision_metrics.putative_int_pyr_pairs = ...
    collision_metrics.rate_difference > rate_change_threshold & preInt_select;

collision_metrics.uLEDResponses_OutInterval  = uLEDResponses_OutInterval;
collision_metrics.uLEDResponses_InInterval   = uLEDResponses_InInterval;

collision_metrics.putative_int_pyr_pairs_list = ...
    [collision_metrics.uLEDResponses_InInterval.presynapticID(collision_metrics.putative_int_pyr_pairs)...
    collision_metrics.uLEDResponses_InInterval.postsynapticID(collision_metrics.putative_int_pyr_pairs)];

pyr = strcmpi(cell_metrics.putativeCellType,'Pyramidal cell');
nw = strcmpi(cell_metrics.putativeCellType,'Narrow Interneuron');
ww = strcmpi(cell_metrics.putativeCellType,'Wide Interneuron');

collision_metrics.inhibitory_connection_probability = size(collision_metrics.putative_int_pyr_pairs_list,1)/(length(find(pyr)) * (length(find(nw)) + length(find(ww))));
collision_metrics.excitatory_connection_probability = size(cell_metrics.putativeConnections.excitatory,1)/(length(find(pyr)) * (length(find(nw)) + length(find(ww))));
collision_metrics.inhibitory_connectionsOut = histcounts(collision_metrics.putative_int_pyr_pairs_list(:,1),[0:length(cell_metrics.UID)]+.5);
collision_metrics.inhibitory_connectionsIn = histcounts(collision_metrics.putative_int_pyr_pairs_list(:,2),[0:length(cell_metrics.UID)]+.5);
collision_metrics.presynapticID = collision_metrics.uLEDResponses_InInterval.presynapticID;
collision_metrics.postsynapticID = collision_metrics.uLEDResponses_InInterval.postsynapticID;

if rand_analysis
    collision_metrics.uLEDResponses_RandInterval   = uLEDResponses_RandInterval;
end 
% % boostraping
for ii = 1:length(collision_metrics.rate_difference)
    rate_difference_rand(ii,:) = uLEDResponses_OutInterval.maxRatePulse(ii) - uLEDResponses_RandInterval.maxRatePulse(ii,:);
    collision_metrics.boostrap_CI_05(ii,:) = prctile(rate_difference_rand(ii,:),[5 97.5]);
    collision_metrics.boostrap_CI_01(ii,:) = prctile(rate_difference_rand(ii,:),[0.5 99.5]);
    collision_metrics.boostrap_CI_001(ii,:) = prctile(rate_difference_rand(ii,:),[0.05 99.95]);
    collision_metrics.boostrap_CI_0001(ii,:) = prctile(rate_difference_rand(ii,:),[0.005 99.995]);
    collision_metrics.boostrap_CI_00001(ii,:) = prctile(rate_difference_rand(ii,:),[0.0005 99.9995]);

    % make boostrap


    collision_metrics.boostrap_CI_05_test(ii) = ~InIntervals(collision_metrics.rate_difference(ii),[collision_metrics.boostrap_CI_05(ii,1)  collision_metrics.boostrap_CI_05(ii,2)]);
    collision_metrics.boostrap_CI_01_test(ii) = ~InIntervals(collision_metrics.rate_difference(ii),[collision_metrics.boostrap_CI_01(ii,1)  collision_metrics.boostrap_CI_01(ii,2)]);
    collision_metrics.boostrap_CI_001_test(ii) = ~InIntervals(collision_metrics.rate_difference(ii),[collision_metrics.boostrap_CI_001(ii,1)  collision_metrics.boostrap_CI_001(ii,2)]);
    collision_metrics.boostrap_CI_0001_test(ii) = ~InIntervals(collision_metrics.rate_difference(ii),[collision_metrics.boostrap_CI_0001(ii,1)  collision_metrics.boostrap_CI_0001(ii,2)]);
end

% select pairs
prePyr_select = collision_metrics.candidate_pyr_pyr_pairs;
preInt_select = collision_metrics.candidate_int_pyr_pairs;

if saveMat
    disp(' Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename, '.', save_as, '.cellinfo.mat'],'collision_metrics','-v7.3');
end

if update_cell_metrics
    cell_metrics.putativeConnections.inhibitory = collision_metrics.putative_int_pyr_pairs_list;
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.cell_metrics.cellinfo.mat'],'cell_metrics');
end


%% plots
if doPlot
    color_pyr = [.9 .3 .3];
    color_int = [.3 .3 .9];
    color_wint = [.3 .9 .9];
    
    prePyr_select = collision_metrics.candidate_pyr_pyr_pairs;
    preInt_select = collision_metrics.candidate_int_pyr_pairs;
    
    try
        figure
        set(gcf,'Position',[200 -500 700 900]);
        nexttile
        hold on
        h1 = histogram(collision_metrics.rate_difference(prePyr_select),...
            [-10:.2:10],'FaceColor',color_pyr,'EdgeColor','none','Normalization','probability');
        h2 = histogram(collision_metrics.rate_difference(preInt_select),...
            [-10:.2:10]+.05,'FaceColor',color_int,'EdgeColor','none','Normalization','probability');
        xlabel('Rate difference (Hz)'); ylabel('Probability');
        set(gca,'TickDir','out');
        legend([h1 h2],'Pre Pyr', 'Pre Int');
        title(strcat('Spike light collision ',{' '}, label));
    catch
        warning('not enough pulses.')   
    end 

    nexttile
    try
        groupCorr(collision_metrics.rate_difference(prePyr_select), collision_metrics.rate_only_light(prePyr_select),...
            'inAxis',true,'MarkerColor',color_pyr);
        groupCorr(collision_metrics.rate_difference(preInt_select), collision_metrics.rate_only_light(preInt_select),...
            'inAxis',true,'MarkerColor',color_int,'labelOffset',2);
        xlabel('Rate difference (Hz)'); ylabel('Rate during light (Hz)');
   catch
        warning('not enough pulses.')   
    end

    try
        nexttile
        hold on
        bar([1 2], [1 1],'FaceColor',[.9 .9 .9],'EdgeColor','none');
        frac = [length(find(collision_metrics.rate_difference(prePyr_select)>rate_change_threshold))/length(find(prePyr_select))];
        bar([1], frac,'FaceColor', color_pyr,'EdgeColor','none');
        text(1,frac+.05,[num2str(round(frac*100,1)), '%'],'Color',color_pyr,'FontSize',12);
        frac = [length(find(collision_metrics.rate_difference(preInt_select)>rate_change_threshold))/length(find(preInt_select))];
        bar([2], frac,'FaceColor', color_int,'EdgeColor','none');
        text(2,frac+.05,[num2str(round(frac*100,1)), '%'],'Color',color_int,'FontSize',12);
        set(gca,'TickDir','out','XTick',[1 2],'XTickLabel',{'Pre pyr','Pre int'},'XTickLabelRotation',45);
    ylabel('Pairs with rate difference > 3 Hz');
   catch
        warning('not enough pulses.')   
    end
    try
        nexttile
        groupStats({(collision_metrics.rateZ_ligh_spike_collision(prePyr_select)),...
            (collision_metrics.rateZ_only_light(prePyr_select)),...
            (collision_metrics.rateZ_ligh_spike_collision(preInt_select)),...
            (collision_metrics.rateZ_only_light(preInt_select))},...
            [1 1 2 2; 1 2 1 2],'plotType','roundPlot','plotData',true,'color',...
            [color_pyr; [.7 .2 .2]; color_int; [.2 .2 .7];],'inAxis',true);
        set(gca,'XTick',[1 2 3.5 4.5],'XTickLabel', {'PrePyr L+Spk', 'PrePyr L', 'PreInt L+Spk', 'PreInt L'},...
            'XTickLabelRotation',45);
        ylabel('Rate responses (Std)');
       catch
        warning('not enough pulses.')   
    end
    try
        nexttile
        groupStats({(collision_metrics.rate_ligh_spike_collision(prePyr_select)),...
            (collision_metrics.rate_only_light(prePyr_select)),...
            (collision_metrics.rate_ligh_spike_collision(preInt_select)),...
            (collision_metrics.rate_only_light(preInt_select))},...
            [1 1 2 2; 1 2 1 2],'plotType','roundPlot','plotData',true,'color',...
            [color_pyr; [.7 .2 .2]; color_int; [.2 .2 .7];],'inAxis',true);
        set(gca,'XTick',[1 2 3.5 4.5],'XTickLabel', {'PrePyr L+Spk', 'PrePyr L', 'PreInt L+Spk', 'PreInt L'},...
            'XTickLabelRotation',45);
        ylabel('Rate responses (Hz)');
       catch
        warning('not enough pulses.')   
    end
    try
        mkdir('SummaryFigures/');
        saveas(gcf,strcat('SummaryFigures\Light_spike_Collision_',label,'.png'));
    catch
        warning('not enough pulses.')   
    end

    % visualize pairs
    try
        figure
        subplot(2,3,[1 2]);
        groupCorr(log10(collision_metrics.rate_only_light(preInt_select)),(collision_metrics.rate_difference(preInt_select)),...
            'inAxis',true,'MarkerColor',color_int,'MarkerSize',15);
        ax = axis;
        plot(ax(1:2), [rate_change_threshold rate_change_threshold], '-r');
        putative_inh_pairs = collision_metrics.putative_int_pyr_pairs;
        plot(log10(collision_metrics.rate_only_light(putative_inh_pairs)), collision_metrics.rate_difference(putative_inh_pairs),'ok');
        xlabel('Rate during light (Hz)'); ylabel('Rate difference [Hz]');
        LogScale('x',10);
        ax1 = ax;
        subplot(2,3,[3]);
        hold on
        [counts, edges] = histcounts(collision_metrics.rate_difference(preInt_select),[-5:.5:15]);
        bins = edges(1:end-1)+ diff(edges(1:2))/2;
        barh(bins, counts,'FaceColor',color_int,'EdgeColor','none','BarWidth',1);
        barh(bins(bins>rate_change_threshold), counts(bins>rate_change_threshold),...
            'FaceColor',[.5 .5 .5],'EdgeColor','none','BarWidth',1,'FaceAlpha',.3);
        ax = axis;
        plot(ax(1:2), [rate_change_threshold rate_change_threshold], '-r');
        ylim(ax1(3:4));
        xlabel('Counts');
        non_putative_inh_pairs = collision_metrics.is_lightResponsive & ~putative_inh_pairs;
        subplot(2,2,3);
        hold on
        plotFill(collision_metrics.uLEDResponses_InInterval.timestamps,...
        collision_metrics.uLEDResponses_InInterval.responsecurve(putative_inh_pairs,:),'color',[.1 .1 .1],'smoothOpt',5,'style','filled');
        plotFill(collision_metrics.uLEDResponses_InInterval.timestamps,...
        collision_metrics.uLEDResponses_OutInterval.responsecurve(putative_inh_pairs,:),'color',[.8 .1 .1],'smoothOpt',5,'style','filled');
        plotFill(collision_metrics.uLEDResponses_InInterval.timestamps,...
        collision_metrics.uLEDResponses_OutInterval.responsecurve(putative_inh_pairs,:) ...
            - collision_metrics.uLEDResponses_InInterval.responsecurve(putative_inh_pairs,:)...
            ,'color',[.1 .1 .1],'smoothOpt',5,'style','edge');
        plotFill(collision_metrics.uLEDResponses_InInterval.timestamps,...
        collision_metrics.uLEDResponses_OutInterval.responsecurve(non_putative_inh_pairs,:) ...
            - collision_metrics.uLEDResponses_InInterval.responsecurve(non_putative_inh_pairs,:)...
            ,'color',[.9 .9 .1],'smoothOpt',5,'style','edge');
        xlim([-.02 .04]);
        ylabel('Hz'); xlabel('Time (s)');
        subplot(2,2,4);
        hold on
        scatter(cell_metrics.troughToPeak(pyr),cell_metrics.burstIndex_Royer2012(pyr),40,'filled',...
            'MarkerFaceColor',color_pyr,'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
        scatter(cell_metrics.troughToPeak(nw),cell_metrics.burstIndex_Royer2012(nw),40,'filled',...
            'MarkerFaceColor',color_int,'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
        scatter(cell_metrics.troughToPeak(ww),cell_metrics.burstIndex_Royer2012(ww),40,'filled',...
            'MarkerFaceColor',color_wint,'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    
        for ii = 1:length(cell_metrics.putativeConnections.excitatory)
            pre_post = [cell_metrics.putativeConnections.excitatory(ii,1) cell_metrics.putativeConnections.excitatory(ii,2)];
            plot(cell_metrics.troughToPeak(pre_post),...
                cell_metrics.burstIndex_Royer2012(pre_post),'-','color',[.4 .1 .1 .2]);
        end
    
        for ii = 1:size(collision_metrics.putative_int_pyr_pairs_list,1)
            pre_post = [collision_metrics.putative_int_pyr_pairs_list(ii,1) collision_metrics.putative_int_pyr_pairs_list(ii,2)];
            plot(cell_metrics.troughToPeak(pre_post),...
                cell_metrics.burstIndex_Royer2012(pre_post),'-','color',[.1 .1 .4]);
        end
        text(0.1, .9,[num2str(size(cell_metrics.putativeConnections.excitatory,1)) 'exc pairs (' num2str(round(collision_metrics.excitatory_connection_probability,2)) ' prob)'],"Units","normalized");
        text(0.1, .8,[num2str(size(collision_metrics.putative_int_pyr_pairs_list,1)) 'inh pairs (' num2str(round(collision_metrics.inhibitory_connection_probability,2)) ' prob)'],"Units","normalized");
        ylabel('Bursting Index (Royer2012)'); xlabel('Through to peak (ms)');
        
        mkdir('SummaryFigures/');
        saveas(gcf,strcat('SummaryFigures\Light_spike_Collision_pairs_',label,'.png'));

    catch
        warning('not enough pulses.')   
    end    
    
    % plotting pairs
    try 
        figure,
        subplot(2,3,1)
        imagesc_ranked(uLEDResponses_OutInterval.timestamps, [], uLEDResponses_OutInterval.responsecurveZ(putative_inh_pairs,:),[-30 30],collision_metrics.rateZ_difference(putative_inh_pairs),'smoothOpt',5);
        ylabel('Z-scored pairs (S.D.)');
        subplot(2,3,2)
        imagesc_ranked(uLEDResponses_OutInterval.timestamps, [], uLEDResponses_InInterval.responsecurveZ(putative_inh_pairs,:),[-30 30],collision_metrics.rateZ_difference(putative_inh_pairs),'smoothOpt',5);
        subplot(2,3,3)
        imagesc_ranked(uLEDResponses_InInterval.timestamps, [], uLEDResponses_OutInterval.responsecurveZ(putative_inh_pairs,:) - uLEDResponses_InInterval.responsecurveZ(putative_inh_pairs,:),...
            [-30 30],collision_metrics.rateZ_difference(putative_inh_pairs),'smoothOpt',5);
    
        subplot(2,3,4)
        imagesc_ranked(uLEDResponses_OutInterval.timestamps, [], uLEDResponses_OutInterval.responsecurve(putative_inh_pairs,:),[0 30],collision_metrics.rate_difference(putative_inh_pairs),'smoothOpt',5);
        ylabel('Rate pairs (Hz)');
        subplot(2,3,5)
        imagesc_ranked(uLEDResponses_OutInterval.timestamps, [], uLEDResponses_InInterval.responsecurve(putative_inh_pairs,:),[0 30],collision_metrics.rate_difference(putative_inh_pairs),'smoothOpt',5);
        xlabel('Time since light stimulation (ms)')
        subplot(2,3,6)
        imagesc_ranked(uLEDResponses_InInterval.timestamps, [], uLEDResponses_OutInterval.responsecurve(putative_inh_pairs,:) - uLEDResponses_InInterval.responsecurve(putative_inh_pairs,:),...
            [0 30],collision_metrics.rate_difference(putative_inh_pairs),'smoothOpt',5);
        colormap jet
        mkdir('SummaryFigures/');
        saveas(gcf,strcat('SummaryFigures\Light_spike_Collision_ranked_pairs_',label,'.png'));
     catch
        warning('not enough pulses.')   
    end
end

cd(prevPath);
end
