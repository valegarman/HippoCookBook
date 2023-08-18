
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
addParameter(p,'rate_change_threshold',5,@isnumeric); % 
addParameter(p,'spikes',[],@isstruct); % 
addParameter(p,'interpolate_pulse_sides',true,@islogical); % 

parse(p, uLEDResponses_interval, varargin{:});
basepath = p.Results.basepath;
cell_metrics = p.Results.cell_metrics;
light_response_threshold_factor = p.Results.light_response_threshold_factor;
light_response_only_pyr = p.Results.light_response_only_pyr;
before_pulse_factor = p.Results.before_pulse_factor;
doPlot = p.Results.doPlot;
label = p.Results.label;
saveMat = p.Results.saveMat;
rate_change_threshold = p.Results.rate_change_threshold;
spikes = p.Results.spikes;
interpolate_pulse_sides = p.Results.interpolate_pulse_sides;

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

uLEDResponses_OutInterval.presynapticID = [];
uLEDResponses_OutInterval.presynapticCellType = [];
uLEDResponses_OutInterval.postsynapticID = [];
uLEDResponses_OutInterval.maxRateBeforePulse = [];
uLEDResponses_OutInterval.maxRatePulse = [];
uLEDResponses_OutInterval.maxZBeforePulse = [];
uLEDResponses_OutInterval.maxZPulse = [];
uLEDResponses_OutInterval.putativeCellType = [];
uLEDResponses_OutInterval.is_in_rateBeforePuse_similar = [];

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
end

% including responseCurves
in_responsecurve = [];
out_responsecurve = [];
in_responsecurveZ = [];
out_responsecurveZ = [];
timestamps = uLEDResponses_interval{1}.in_interval.timestamps;
for ii = 1:length(uLEDResponses_interval)
    for jj = 1:length(uLEDResponses_interval{ii}.in_interval.maxRespLED.LEDs)
        maxLEDs = uLEDResponses_interval{ii}.in_interval.maxRespLED.LEDs(jj);
        if isnan(maxLEDs) || ii == jj
            in_responsecurve = [in_responsecurve; timestamps'*NaN];
            out_responsecurve = [out_responsecurve; timestamps'*NaN];
            in_responsecurveZ = [in_responsecurveZ; timestamps'*NaN];
            out_responsecurveZ = [out_responsecurveZ; timestamps'*NaN];
        else
            in_responsecurve = ...
                [in_responsecurve; squeeze(uLEDResponses_interval{ii}.in_interval.responsecurve(jj,1,maxLEDs,:))'];
            out_responsecurve = ...
                [out_responsecurve; squeeze(uLEDResponses_interval{ii}.out_interval.responsecurve(jj,1,maxLEDs,:))'];
            in_responsecurveZ = ...
                [in_responsecurveZ; squeeze(uLEDResponses_interval{ii}.in_interval.responsecurveZ(jj,1,maxLEDs,:))'];
            out_responsecurveZ = ...
                [out_responsecurveZ; squeeze(uLEDResponses_interval{ii}.out_interval.responsecurveZ(jj,1,maxLEDs,:))'];
        end
    end
end
uLEDResponses_OutInterval.responsecurve = out_responsecurve;
uLEDResponses_OutInterval.responsecurveZ = out_responsecurveZ;
uLEDResponses_OutInterval.timestamps = timestamps';
uLEDResponses_InInterval.responsecurve = in_responsecurve;
uLEDResponses_InInterval.responsecurveZ = in_responsecurveZ;
uLEDResponses_InInterval.timestamps = timestamps';

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
            uLEDResponses_OutInterval.responsecurveZ(jj,samples_to_interpolate(ii,:)) = ...
                interp1(x_axis,uLEDResponses_OutInterval.responsecurveZ(jj,x_axis),samples_to_interpolate(ii,:));
            uLEDResponses_InInterval.responsecurve(jj,samples_to_interpolate(ii,:)) = ...
                interp1(x_axis,uLEDResponses_InInterval.responsecurve(jj,x_axis),samples_to_interpolate(ii,:));
            uLEDResponses_InInterval.responsecurveZ(jj,samples_to_interpolate(ii,:)) = ...
                interp1(x_axis,uLEDResponses_InInterval.responsecurveZ(jj,x_axis),samples_to_interpolate(ii,:));
        end
    end
end

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
        collision_metrics.is_lightResponsive & collision_metrics.is_rate_before_similar_factor & ...
        collision_metrics.is_post_pyramidalCell & not_same_cell;
collision_metrics.candidate_pyr_pyr_pairs = collision_metrics.is_pre_interneuron==0 & ...
        collision_metrics.is_lightResponsive & collision_metrics.is_rate_before_similar_factor & ...
        collision_metrics.is_post_pyramidalCell & not_same_cell;

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
collision_metrics.pre_waveforms   = pre_waveforms;
collision_metrics.post_waveforms  = post_waveforms;
collision_metrics.pre_firingRate  = pre_firingRate;
collision_metrics.post_firingRate = post_firingRate;

% actual metrics
collision_metrics.rate_difference            = uLEDResponses_OutInterval.maxRatePulse - uLEDResponses_InInterval.maxRatePulse;
collision_metrics.rateZ_difference           = uLEDResponses_OutInterval.maxZPulse - uLEDResponses_InInterval.maxZPulse;
collision_metrics.rate_only_light            = uLEDResponses_OutInterval.maxRatePulse;
collision_metrics.rateZ_only_light           = uLEDResponses_OutInterval.maxZPulse;
collision_metrics.rate_ligh_spike_collision  = uLEDResponses_InInterval.maxRatePulse;
collision_metrics.rateZ_ligh_spike_collision = uLEDResponses_InInterval.maxZPulse;
% 

collision_metrics.uLEDResponses_OutInterval  = uLEDResponses_OutInterval;
collision_metrics.uLEDResponses_InInterval   = uLEDResponses_InInterval;

% select pairs
figure
nexttile
groupCorr(collision_metrics.rate_only_light(prePyr_select),collision_metrics.rate_difference(prePyr_select),...
        'inAxis',true,'MarkerColor',color_pyr,'labelOffset',2);
groupCorr(collision_metrics.rate_only_light(preInt_select),collision_metrics.rate_difference(preInt_select),...
        'inAxis',true,'MarkerColor',color_int,'labelOffset',2);
ylabel('Rate difference (Hz)'); xlabel('Rate during light (Hz)');

nexttile
groupCorr(collision_metrics.rateZ_only_light(prePyr_select),collision_metrics.rateZ_difference(prePyr_select),...
        'inAxis',true,'MarkerColor',color_pyr,'labelOffset',2);
groupCorr(collision_metrics.rateZ_only_light(preInt_select),collision_metrics.rateZ_difference(preInt_select),...
        'inAxis',true,'MarkerColor',color_int,'labelOffset',2);
ylabel('Rate difference (Std)'); xlabel('Rate during light (Std)');

collision_metrics.rate_difference

if saveMat
    disp(' Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.lightSpikeCollisions.cellinfo.mat'],'collision_metrics','-v7.3');
end

%% plots
if doPlot
    color_pyr = [.9 .3 .3];
    color_int = [.3 .3 .9];
    
    prePyr_select = collision_metrics.candidate_pyr_pyr_pairs;
    preInt_select = collision_metrics.candidate_int_pyr_pairs;
    
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
    
    nexttile
    groupCorr(collision_metrics.rate_difference(prePyr_select), collision_metrics.rate_only_light(prePyr_select),...
        'inAxis',true,'MarkerColor',color_pyr);
    groupCorr(collision_metrics.rate_difference(preInt_select), collision_metrics.rate_only_light(preInt_select),...
        'inAxis',true,'MarkerColor',color_int,'labelOffset',2);
    xlabel('Rate difference (Hz)'); ylabel('Rate during light (Hz)');

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
    
    mkdir('SummaryFigures/')
    saveas(gcf,strcat('SummaryFigures\Light_spike_Collision_',label,'.png'));

    % visualize pairs
    
    

    
end

cd(prevPath);
end