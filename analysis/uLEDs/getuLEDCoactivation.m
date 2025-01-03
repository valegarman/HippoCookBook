
function uLEDcoactivation = getuLEDcoactivation(varargin)
% Compute several measures for uled coactivation
%
% To do: include distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'uledpulses',getuLEDPulses,@isstruct);
addParameter(p,'session',loadSession,@isstruct);
addParameter(p,'cell_metrics',loadCellMetrics,@isstruct);
addParameter(p,'spikes',loadSpikes,@isstruct);
addParameter(p,'winSize_uled',.2,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'save_as','uLEDcoactivation',@ischar);
addParameter(p,'before_pulse_win',[-.1 -0.05],@isnumeric);
addParameter(p,'during_pulse_win',[0.01 .020],@isnumeric);
addParameter(p,'winSizePlot',[-.1 .1],@isnumeric);
addParameter(p,'winCoactivation',.01,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
uledpulses = p.Results.uledpulses;
session = p.Results.session;
cell_metrics = p.Results.cell_metrics;
spikes = p.Results.spikes;
winSize_uled = p.Results.winSize_uled;
doPlot = p.Results.doPlot;
force = p.Results.force;
saveMat = p.Results.saveMat;
save_as = p.Results.save_as;
before_pulse_win = p.Results.before_pulse_win;
during_pulse_win = p.Results.during_pulse_win;
winSizePlot = p.Results.winSizePlot;
winCoactivation = p.Results.winCoactivation;

prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDcoactivation.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('uLED coactivation already computed! Loading file...');
    load(targetFile.name);
    return
end

%% COACTIVATION ANALYSIS
% 0 Compute uled responses according to local parameters
uledResponses = getuLEDResponse('force', true, 'winSize', winSize_uled, 'before_pulse_win',before_pulse_win, 'during_pulse_win', during_pulse_win, 'winSizePlot', winSizePlot, 'saveMat', false, 'save_as', 'uledResponses_coactivation');
uLEDcoactivation.uledResponses = uledResponses;

% 1.1 Compute CCG 
epochsNames = [];
epochsInts = [];
for ii = 1:size(session.epochs,2)
    epochsNames{ii} = session.epochs{ii}.behavioralParadigm;
    epochsInts(ii,:) = [session.epochs{ii}.startTime session.epochs{ii}.stopTime];
end
if sum(ismember({'precoactivation', 'coactivation', 'postcoactivation'}, epochsNames)) ~= 3
    error('Epochs names do not match the expected results! Define "precoactivation", "coactivation", and "postcoactivation" epochs in the current session!');
end

disp('Computing CCG...');
disp('... CCG during pre-coactivation');
averageCCG_precoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'precoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during coactivation');
averageCCG_coact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'coactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during post-coactivation');
averageCCG_postcoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'postcoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);

uLEDcoactivation.epochsNames = epochsNames;
uLEDcoactivation.epochsInts = epochsInts;
uLEDcoactivation.averageCCG_precoact = averageCCG_precoact;
uLEDcoactivation.averageCCG_coact = averageCCG_coact;
uLEDcoactivation.averageCCG_postcoact = averageCCG_postcoact;

% distance 
coordinates = [cell_metrics.trilat_x', cell_metrics.trilat_y'];
D = pdist2(coordinates, coordinates);
distance.units = D;
distance.pairs = reshape(D',[],1);
clear D
maxWaveformCh1 = cell_metrics.maxWaveformCh1';
D = pdist2(maxWaveformCh1, maxWaveformCh1);
same_electrode = D == 0;
distance.same_electrode_units = same_electrode;
distance.same_electrode_pairs = reshape(same_electrode',[],1);

% Compute z-score coactivation before, during and after
units.preZ = [];
units.coactZ = [];
units.postZ = [];
pairs.curve_preZ = [];
pairs.curve_coactZ = [];
pairs.curve_postZ = [];
pairs.pre_id = [];
pairs.post_id = [];

zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
winCoactivation = InIntervals(averageCCG_precoact.timestamps, [-winCoactivation winCoactivation]);
win_Z = InIntervals(averageCCG_precoact.timestamps, [before_pulse_win(1) before_pulse_win(2)]);
for kk = 1:size(averageCCG_precoact.allCcg,2)
    for jj = 1:size(averageCCG_precoact.allCcg,2)
        pairs.pre_id = [pairs.pre_id; kk];
        pairs.post_id = [pairs.post_id; jj];
        if kk == jj
            units.preZ(kk,jj) = NaN;
            units.coactZ(kk,jj) = NaN;
            units.postZ(kk,jj) = NaN;
            
            temp = nan(size(averageCCG_precoact.allCcg(:,1,1)));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            pairs.curve_postZ = [pairs.curve_postZ; temp'];
        else
            % units
            % pre
            temp = averageCCG_precoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.preZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            % coact
            temp = averageCCG_coact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.coactZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            % post
            temp = averageCCG_postcoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.postZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_postZ = [pairs.curve_postZ; temp'];
        end
    end
end
winCoactivation(zero_ind) = 0;
pairs.preZ = nanmean(pairs.curve_preZ(:,winCoactivation),2);
pairs.coactZ = nanmean(pairs.curve_coactZ(:,winCoactivation),2);
pairs.postZ = nanmean(pairs.curve_postZ(:,winCoactivation),2);

% remove inf
pairs.preZ(find(isinf(pairs.preZ))) = NaN;
pairs.coactZ(find(isinf(pairs.coactZ))) = NaN;
pairs.postZ(find(isinf(pairs.postZ))) = NaN;

units.preZ(find(isinf(units.preZ))) = NaN;
units.coactZ(find(isinf(units.coactZ))) = NaN;
units.postZ(find(isinf(units.postZ))) = NaN;

uLEDcoactivation.units = units;
uLEDcoactivation.pairs = pairs;

% pyramidal neurons (label as NaN all non pyramidal cells)
pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
uLEDcoactivation.units_pyramidalCells = label_unit_as_nan(uLEDcoactivation.units, ~pyramidal_neurons);

pyramidal_neuron_pair = ismember(uLEDcoactivation.pairs.pre_id, find(pyramidal_neurons)) & ismember(uLEDcoactivation.pairs.post_id, find(pyramidal_neurons));
uLEDcoactivation.pairs_pyramidalCells = label_pairs_as_nan(uLEDcoactivation.pairs, ~pyramidal_neuron_pair);

% light responsive neurons and pyramidal neurons
boostrap_mat = nansum(squeeze(abs(uledResponses.bootsTrapTest(:,1,:))),2)>0;
targetCells = pyramidal_neurons & boostrap_mat;
uLEDcoactivation.units_pyramidalCells_lightResponsive = label_unit_as_nan(uLEDcoactivation.units, ~targetCells);

lightResponsive_pyr_pair = ismember(uLEDcoactivation.pairs.pre_id, find(targetCells)) & ismember(uLEDcoactivation.pairs.post_id, find(targetCells));
uLEDcoactivation.pairs_pyramidalCells_lightResponsive = label_pairs_as_nan(uLEDcoactivation.pairs, ~lightResponsive_pyr_pair);

% light responsive neurons
boostrap_mat = nansum(squeeze(abs(uledResponses.bootsTrapTest(:,1,:))),2)>0;
targetCells = boostrap_mat;
uLEDcoactivation.units_lightResponsive = label_unit_as_nan(uLEDcoactivation.units, ~targetCells);

lightResponsive_pair = ismember(uLEDcoactivation.pairs.pre_id, find(targetCells)) & ismember(uLEDcoactivation.pairs.post_id, find(targetCells));
uLEDcoactivation.pairs_lightResponsive = label_pairs_as_nan(uLEDcoactivation.pairs, ~lightResponsive_pair);

% Neurons in different electrode
uLEDcoactivation.units_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.units, distance.same_electrode_units);
uLEDcoactivation.pair_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs, distance.same_electrode_pairs);

% Pyramidal neurons in different electrode
uLEDcoactivation.units_pyramidalCells_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.units_pyramidalCells, distance.same_electrode_units);
uLEDcoactivation.pair_pyramidalCells_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells, distance.same_electrode_pairs);

% Light responsive, pyramidal neurons, in different electrode
uLEDcoactivation.units_pyramidalCells_lightResponsive_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.units_pyramidalCells_lightResponsive, distance.same_electrode_units);
uLEDcoactivation.pairs_pyramidalCells_lightResponsive_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells_lightResponsive, distance.same_electrode_pairs);

% Light responsive
uLEDcoactivation.units_lightResponsive_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.units_lightResponsive, distance.same_electrode_units);
uLEDcoactivation.pairs_lightResponsive_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs_lightResponsive, distance.same_electrode_pairs);

% 1.2 % Coactivated LEDs
% get pairs of coactivated uleds
[status] = InIntervals(uledpulses.timestamps(:,1),epochsInts(find(ismember(epochsNames,'coactivation')),:));
times_coactivation = uledpulses.timestamps(status==1,1);
codes_coactivation = uledpulses.code(status==1,1);
% list_codes = unique(codes_coactivation);
uLEDcoactivation.coactivation_times = [];
for kk = 1:12
    uLEDcoactivation.coactivation_times{kk} = times_coactivation(find(codes_coactivation==kk));
end
binSize = [0.001];
winSize = [1];
[allCcg, t_ccg] = CCG(uLEDcoactivation.coactivation_times,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);
winCoactivation = p.Results.winCoactivation;
winCoactivation = InIntervals(t_ccg, [-winCoactivation winCoactivation]);
for kk = 1:size(allCcg,2)
    for jj = 1:size(allCcg,2)
        temp = zscore(allCcg(:,kk,jj));
        coactivation_matrix(kk,jj) = max(temp(winCoactivation)) > 10;
    end
end
uLEDcoactivation.uled_coactivation_matrix = coactivation_matrix;

% classifying neurons according to response
zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
winCoactivation = p.Results.winCoactivation;
winCoactivation = InIntervals(averageCCG_precoact.timestamps, [-winCoactivation winCoactivation]);
win_Z = InIntervals(averageCCG_precoact.timestamps, [before_pulse_win(1) before_pulse_win(2)]);
boostrap_mat = squeeze(uledResponses.bootsTrapTest(:,1,:));

% coactivation zscore for uleds using light responsive pyramidal neurons
uleds_temp.preZ = [];
uleds_temp.coactZ = [];
uleds_temp.postZ = [];
uleds_temp.post_pre = [];
uleds_temp.rateZmat = [];
for kk = 1:size(boostrap_mat,2)
    for jj = 1:size(boostrap_mat,2)
        % find neurons
        neurons_kk = find(boostrap_mat(:,kk)==1 & pyramidal_neurons);
        neurons_jj  = find(boostrap_mat(:,jj )==1 & pyramidal_neurons);
        uleds_temp.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

        % coactivations
        temp = units.preZ(neurons_kk,neurons_jj);
        uleds_temp.preZ(kk,jj) = nanmean(temp(:));
        temp = units.coactZ(neurons_kk,neurons_jj);
        uleds_temp.coactZ(kk,jj) = nanmean(temp(:));
        temp = units.postZ(neurons_kk,neurons_jj);
        uleds_temp.postZ(kk,jj) = nanmean(temp(:));
        uleds_temp.post_pre(kk,jj) = uleds_temp.postZ(kk,jj) - uleds_temp.preZ(kk,jj);
    end
end
uLEDcoactivation.uleds_pyramidalCells_lightResponsive = uleds_temp;

% coactivation zscore for uleds using light responsive
uleds_temp.preZ = [];
uleds_temp.coactZ = [];
uleds_temp.postZ = [];
uleds_temp.post_pre = [];
uleds_temp.rateZmat = [];
for kk = 1:size(boostrap_mat,2)
    for jj = 1:size(boostrap_mat,2)
        % find neurons
        neurons_kk = find(boostrap_mat(:,kk)==1);
        neurons_jj  = find(boostrap_mat(:,jj )==1);
        uleds_temp.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

        % coactivations
        temp = units.preZ(neurons_kk,neurons_jj);
        uleds_temp.preZ(kk,jj) = nanmean(temp(:));
        temp = units.coactZ(neurons_kk,neurons_jj);
        uleds_temp.coactZ(kk,jj) = nanmean(temp(:));
        temp = units.postZ(neurons_kk,neurons_jj);
        uleds_temp.postZ(kk,jj) = nanmean(temp(:));
        uleds_temp.post_pre(kk,jj) = uleds_temp.postZ(kk,jj) - uleds_temp.preZ(kk,jj);
    end
end
uLEDcoactivation.uleds_lightResponsive = uleds_temp;


% 1.3 % explained variance
evStats_units = explained_variance(spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~pyramidal_neurons) = [];
evStats_pyramidalCells = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~targetCells) = [];
evStats_pyramidalCells_lightResponsive = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
uLEDcoactivation.evStats_units = evStats_units;
uLEDcoactivation.evStats_pyramidalCells = evStats_pyramidalCells;
uLEDcoactivation.evStats_pyramidalCells_lightResponsive = evStats_pyramidalCells_lightResponsive;

% Saving 
if saveMat
    disp('Saving...');
    save([basenameFromBasepath(pwd) '.' save_as '.cellinfo.mat'],'uLEDcoactivation');
end

if doPlot
    
    % 1 stats
    % all
    figure('Position', [100, 100, 1200, 1000]); % Create and set size in one step
    tiledlayout(3,6)
    nexttile
    groupStats({uLEDcoactivation.pairs.preZ, uLEDcoactivation.pairs.coactZ, uLEDcoactivation.pairs.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylabel('Cofiring (SD)');
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim([-3 3]);
    title(['All neurons'],'FontWeight','normal');
    
    % Pyramidal neurons
    nexttile
    groupStats({uLEDcoactivation.pairs_pyramidalCells.preZ, uLEDcoactivation.pairs_pyramidalCells.coactZ, uLEDcoactivation.pairs_pyramidalCells.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim([-3 3]);
    title(['Pyr neurons'],'FontWeight','normal');

    % Light responsive pyr
    nexttile
    groupStats({uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ, uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ, uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim([-3 3]);
    title(['Light responsive Pyr'],'FontWeight','normal');

    % Abs all
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs.preZ)), log10(abs(uLEDcoactivation.pairs.coactZ)), log10(abs(uLEDcoactivation.pairs.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim([-2 2]);
    LogScale('y',10);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    title(['All neurons'],'FontWeight','normal');
    ylabel('abs cofiring (SD)');

    % Abs all
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs_pyramidalCells.preZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim([-2 2]);
    LogScale('y',10);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    title(['Pyr neurons'],'FontWeight','normal');
    ylabel('abs cofiring (SD)');

    % Abs pyr
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim([-2 2]);
    LogScale('y',10);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    title(['Light responsive Pyr'],'FontWeight','normal');
    ylabel('abs cofiring (SD)');

    % 2 corr

    % correlation all neurons
    nexttile
    x = sign(uLEDcoactivation.pairs.coactZ) .* log10(abs(uLEDcoactivation.pairs.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) .* log10(abs(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coactivation cofiring (log10 Z)');
    ylabel('Abs post-pre (log10 Z)');

    % correlation pyramidal neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % correlation light responsive pyramidal neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

     % correlation abs all neurons
    nexttile
    x = sign(uLEDcoactivation.pairs.coactZ) .* log10(abs(uLEDcoactivation.pairs.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) .* log10(abs(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coactivation cofiring (log10 Z)');
    ylabel('Abs post-pre (log10 Z)');

    % correlation abs pyramidal neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % correlation abs light responsive pyramidal neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % 2 positive and negative corr
    nexttile
    x = sign(uLEDcoactivation.pairs.coactZ) .* log10(abs(uLEDcoactivation.pairs.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) .* log10(abs(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    % 3 explained variance
    nexttile
    imagesc(1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), 1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), uLEDcoactivation.evStats_pyramidalCells.P1, [-.5 .5]);
    colormap jet
    title(['Pre'],'FontWeight','normal');
    axis square
    ylabel('Pyr neurons');
    xlabel('Pyr neurons');
        
    nexttile
    imagesc(1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), 1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), uLEDcoactivation.evStats_pyramidalCells.R, [-.5 .5]);
    colormap jet
    title(['Coact'],'FontWeight','normal');
    axis square
    xlabel({['Explain variance (EV vs REV)']; ['all: ' num2str(round(uLEDcoactivation.evStats_units.ev,1)) ' vs ' num2str(round(uLEDcoactivation.evStats_units.rev,1))]; ...
        ['Pyr: ' num2str(round(uLEDcoactivation.evStats_pyramidalCells.ev,1)) ' vs ' num2str(round(uLEDcoactivation.evStats_pyramidalCells.rev,1))]; ...
        ['Light resp pyr: ' num2str(round(uLEDcoactivation.evStats_pyramidalCells_lightResponsive.ev,1)) ' vs ' num2str(round(uLEDcoactivation.evStats_pyramidalCells_lightResponsive.rev,1))]});

    nexttile
    imagesc(1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), 1:size(uLEDcoactivation.evStats_pyramidalCells.P1,1), uLEDcoactivation.evStats_pyramidalCells.P2, [-.5 .5]);
    colormap jet
    title(['Post'],'FontWeight','normal');
    axis square
    
    exportgraphics(gcf,['SummaryFigures\', save_as, '.png']);
end

cd(prevPath);

end

function output_struct = label_pairs_as_nan(input_struct, who_is_nan)
% input_sctruct: structure with some with data that should be
% labeled as NaN
output_struct = input_struct;
fnames = fieldnames(output_struct);
for ii = 1:length(fnames)
    output_struct.(fnames{ii})(who_is_nan,:) = NaN;
end

end

function output_struct = label_unit_as_nan(input_struct, who_is_nan)
% input_sctruct: structure with some fields with pairwise data that should be
% labeled as NaN column-wise and row-wise
output_struct = input_struct;
fnames = fieldnames(output_struct);
for ii = 1:length(fnames)
    output_struct.(fnames{ii})(who_is_nan,:) = NaN;
    output_struct.(fnames{ii})(:,who_is_nan) = NaN;
end

end

function output_struct = mask_unit_as_nan(input_struct, nan_mask)
% input_sctruct: structure with some fields with pairwise data that should be
% labeled as NaN from a same size binary mask
output_struct = input_struct;
fnames = fieldnames(output_struct);
for ii = 1:length(fnames)
    output_struct.(fnames{ii})(nan_mask) = NaN;
end

end

