
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
addParameter(p,'before_pulse_win',[-0.1 -0.06],@isnumeric);
addParameter(p,'during_pulse_win',[0.005 0.045],@isnumeric);
addParameter(p,'winSizePlot',[-.1 .1],@isnumeric);
addParameter(p,'winCoactivation',.01,@isnumeric);
addParameter(p,'winSTD',[-.3 -.15],@isnumeric);
addParameter(p,'offset_precoactivation', 0,@isnumeric);
addParameter(p,'test_uled_responses', 'zscore',@ischar);
addParameter(p,'numRep_uledResp', 500, @isscalar);
addParameter(p,'curate_responses', true, @islogical);
% addParameter(p,'epochs_names', {'precoactivation', 'coactivation', 'postcoactivation'} ,@iscell);

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
offset_precoactivation = p.Results.offset_precoactivation;
test_uled_responses = p.Results.test_uled_responses;
winSTD = p.Results.winSTD;
numRep_uledResp = p.Results.numRep_uledResp;
curate_responses = p.Results.curate_responses;
% epochs_names = p.Results.epochs_names;

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
uledResponses = getuLEDResponse('force', true, 'winSize', winSize_uled, 'before_pulse_win',before_pulse_win, 'during_pulse_win', during_pulse_win, 'winSizePlot', winSizePlot, 'saveMat', false,...
    'save_as', 'uledResponses_coactivation', 'test_in_plots',test_uled_responses, 'zscore_threshold', 1, 'numRep', numRep_uledResp);
switch lower(test_uled_responses)
    case 'boostrap'
        lightRespMat = squeeze((uledResponses.bootsTrapTest(:,1,:)));
    case 'zscore'
        lightRespMat = squeeze((uledResponses.zscoreTest(:,1,:)));
    case 'salt'
        lightRespMat = squeeze((uledResponses.salt.p_value(:,1,:)))<0.05;
    case 'kolmogorov'
        lightRespMat = squeeze((uledResponses.modulationSignificanceLevel(:,1,:)))<0.05;
    otherwise
        warning('Test for defining light responsive neurons do not recognized! Using boostrap results...');
        lightRespMat = squeeze((uledResponses.bootsTrapTest(:,1,:)));
end
lightRespNeurons = nansum(lightRespMat==1,2)>0;

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
averageCCG_precoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'precoactivation')),:) + [0 offset_precoactivation], 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during coactivation');
averageCCG_coact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'coactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during post-coactivation');
averageCCG_postcoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'postcoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);

uLEDcoactivation.epochsNames = epochsNames;
uLEDcoactivation.epochsInts = epochsInts;
uLEDcoactivation.averageCCG_precoact = averageCCG_precoact;
uLEDcoactivation.averageCCG_coact = averageCCG_coact;
uLEDcoactivation.averageCCG_postcoact = averageCCG_postcoact;

% distance of neurons
coordinates = [cell_metrics.trilat_x', cell_metrics.trilat_y'];
D = pdist2(coordinates, coordinates);
distance.pairwise = D;
distance.pairs = reshape(D',[],1);
clear D
maxWaveformCh1 = cell_metrics.maxWaveformCh1';
D = pdist2(maxWaveformCh1, maxWaveformCh1);
same_electrode = D == 0;
distance.are_same_electrode_pairwise = same_electrode;
distance.are_same_electrode_pairs = reshape(same_electrode',[],1);
distance.are_distint_electrode_pairwise = distance.are_same_electrode_pairwise == 0;
distance.are_distint_electrode_pairs = distance.are_same_electrode_pairs == 0;

% shanks
same_shank = zeros(spikes.numcells);
for ii = 1:length(spikes.shankID)
    for jj = 1:length(spikes.shankID)
        if spikes.shankID(ii) == spikes.shankID(jj)
            same_shank(ii,jj) = 1;
        end
    end
end
distance.are_same_shank_pairwise = logical(same_shank);
distance.are_same_shank_pairs = logical(reshape(same_shank',[],1));
distance.are_distint_shank_pairwise = distance.are_same_shank_pairwise == 0;
distance.are_distint_shank_pairs = distance.are_same_shank_pairs == 0;

uLEDcoactivation.distance = distance;

% Compute z-score coactivation before, during and after
pairwise.preZ = [];
pairwise.coactZ = [];
pairwise.postZ = [];
pairs.curve_preZ = [];
pairs.curve_coactZ = [];
pairs.curve_postZ = [];
pairs.pre_id = [];
pairs.post_id = [];

zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
winCoactivation = InIntervals(averageCCG_precoact.timestamps, [-winCoactivation winCoactivation]);
win_Z = InIntervals(averageCCG_precoact.timestamps, [winSTD(1) winSTD(2)]);
for kk = 1:size(averageCCG_precoact.allCcg,2)
    for jj = 1:size(averageCCG_precoact.allCcg,2)
        pairs.pre_id = [pairs.pre_id; kk];
        pairs.post_id = [pairs.post_id; jj];
        if kk == jj
            pairwise.preZ(kk,jj) = NaN;
            pairwise.coactZ(kk,jj) = NaN;
            pairwise.postZ(kk,jj) = NaN;
            
            temp = nan(size(averageCCG_precoact.allCcg(:,1,1)));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            pairs.curve_postZ = [pairs.curve_postZ; temp'];
        else
            % pairwise
            % pre
            temp = averageCCG_precoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            pairwise.preZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            % coact
            temp = averageCCG_coact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            pairwise.coactZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            % post
            temp = averageCCG_postcoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            pairwise.postZ(kk,jj) = nanmean(temp(winCoactivation));
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

pairwise.preZ(find(isinf(pairwise.preZ))) = NaN;
pairwise.coactZ(find(isinf(pairwise.coactZ))) = NaN;
pairwise.postZ(find(isinf(pairwise.postZ))) = NaN;

uLEDcoactivation.pairwise = pairwise;
uLEDcoactivation.pairs = pairs;

% pyramidal neurons (label as NaN all non pyramidal cells)
pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
uLEDcoactivation.pairwise_pyramidalCells = label_unit_as_nan(uLEDcoactivation.pairwise, ~pyramidal_neurons);

pyramidal_neuron_pair = ismember(uLEDcoactivation.pairs.pre_id, find(pyramidal_neurons)) & ismember(uLEDcoactivation.pairs.post_id, find(pyramidal_neurons));
uLEDcoactivation.pairs_pyramidalCells = label_pairs_as_nan(uLEDcoactivation.pairs, ~pyramidal_neuron_pair);

uLEDcoactivation.is_pyramidalCells = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';

% light responsive neurons and pyramidal neurons
targetCells = pyramidal_neurons & lightRespNeurons;
uLEDcoactivation.pairwise_pyramidalCells_lightResponsive = label_unit_as_nan(uLEDcoactivation.pairwise, ~targetCells);

lightResponsive_pyr_pair = ismember(uLEDcoactivation.pairs.pre_id, find(targetCells)) & ismember(uLEDcoactivation.pairs.post_id, find(targetCells));
uLEDcoactivation.pairs_pyramidalCells_lightResponsive = label_pairs_as_nan(uLEDcoactivation.pairs, ~lightResponsive_pyr_pair);

uLEDcoactivation.is_pyramidalCells_lightResponsive = pyramidal_neurons & lightRespNeurons;

% light responsive neurons
targetCells = lightRespNeurons;
uLEDcoactivation.pairwise_lightResponsive = label_unit_as_nan(uLEDcoactivation.pairwise, ~targetCells);

lightResponsive_pair = ismember(uLEDcoactivation.pairs.pre_id, find(targetCells)) & ismember(uLEDcoactivation.pairs.post_id, find(targetCells));
uLEDcoactivation.pairs_lightResponsive = label_pairs_as_nan(uLEDcoactivation.pairs, ~lightResponsive_pair);

uLEDcoactivation.is_lightResponsive = lightRespNeurons;

% Neurons in different electrode
uLEDcoactivation.pairwise_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.pairwise, distance.are_same_electrode_pairwise);
uLEDcoactivation.pair_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs, distance.are_same_electrode_pairs);
uLEDcoactivation.is_distinct_electrode = ~distance.are_same_electrode_pairwise;

% Pyramidal neurons in different electrode
uLEDcoactivation.pairwise_pyramidalCells_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells, distance.are_same_electrode_pairwise);
uLEDcoactivation.pair_pyramidalCells_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells, distance.are_same_electrode_pairs);
uLEDcoactivation.is_pyramidalCells_distinct_electrode = pyramidal_neurons & ~distance.are_same_electrode_pairwise;

% Light responsive, pyramidal neurons, in different electrode
uLEDcoactivation.pairwise_pyramidalCells_lightResponsive_distinct_electrode = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_lightResponsive, distance.are_same_electrode_pairwise);
uLEDcoactivation.pairs_pyramidalCells_lightResponsive_distinct_electrode = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells_lightResponsive, distance.are_same_electrode_pairs);
uLEDcoactivation.is_pyramidalCells_lightResponsive_distinct_electrode = pyramidal_neurons & lightRespNeurons & ~distance.are_same_electrode_pairwise;

% pairwise in different shank
uLEDcoactivation.pairwise_distinct_shank = mask_unit_as_nan(uLEDcoactivation.pairwise, distance.are_same_shank_pairwise); % label as Nan those pairwise pairs in same shank
uLEDcoactivation.pair_distinct_shank = label_pairs_as_nan(uLEDcoactivation.pairs, distance.are_same_shank_pairs);
uLEDcoactivation.is_distinct_shank = ~distance.are_same_shank_pairwise;

% Pyramidal neurons in different shank
uLEDcoactivation.pairwise_pyramidalCells_distinct_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells, distance.are_same_shank_pairwise);
uLEDcoactivation.pair_pyramidalCells_distinct_shank = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells, distance.are_same_shank_pairs);
uLEDcoactivation.is_pyramidalCells_distinct_shank = pyramidal_neurons & ~distance.are_same_shank_pairwise;

% Light responsive, pyramidal neurons, in different shank
uLEDcoactivation.pairwise_pyramidalCells_lightResponsive_distinct_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_lightResponsive, distance.are_same_shank_pairwise);
uLEDcoactivation.pairs_pyramidalCells_lightResponsive_distinct_shank = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells_lightResponsive, distance.are_same_shank_pairs);
uLEDcoactivation.is_pyramidalCells_lightResponsive_distinct_shank = pyramidal_neurons & lightRespNeurons & ~distance.are_same_shank_pairwise;

% pairwise in same shank
uLEDcoactivation.pairwise_same_shank = mask_unit_as_nan(uLEDcoactivation.pairwise, distance.are_distint_shank_pairwise); % label as Nan those pairwise pairs in distint shank
uLEDcoactivation.pair_same_shank = label_pairs_as_nan(uLEDcoactivation.pairs, distance.are_distint_shank_pairs);
uLEDcoactivation.is_same_shank = ~distance.are_distint_shank_pairwise;

% Pyramidal neurons in same shank
uLEDcoactivation.pairwise_pyramidalCells_same_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells, distance.are_distint_shank_pairwise);
uLEDcoactivation.pair_pyramidalCells_same_shank = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells, distance.are_distint_shank_pairs);
uLEDcoactivation.is_pyramidalCells_same_shank = pyramidal_neurons & ~distance.are_distint_shank_pairwise;

% Light responsive, pyramidal neurons, in same shank
uLEDcoactivation.pairwise_pyramidalCells_lightResponsive_same_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_lightResponsive, distance.are_distint_shank_pairwise);
uLEDcoactivation.pairs_pyramidalCells_lightResponsive_same_shank = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells_lightResponsive, distance.are_distint_shank_pairs);
uLEDcoactivation.is_pyramidalCells_lightResponsive_same_shank = pyramidal_neurons & lightRespNeurons & ~distance.are_distint_shank_pairwise;

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
% coactivation zscore for uleds using light responsive pyramidal neurons
uleds_temp.preZ = [];
uleds_temp.coactZ = [];
uleds_temp.postZ = [];
uleds_temp.post_pre = [];
uleds_temp.rateZmat = [];
for kk = 1:size(lightRespMat,2)
    for jj = 1:size(lightRespMat,2)
        % find neurons
        neurons_kk = find(lightRespMat(:,kk)==1 & pyramidal_neurons);
        neurons_jj  = find(lightRespMat(:,jj )==1 & pyramidal_neurons);
        uleds_temp.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

        % coactivations
        temp = pairwise.preZ(neurons_kk,neurons_jj);
        uleds_temp.preZ(kk,jj) = nanmean(temp(:));
        temp = pairwise.coactZ(neurons_kk,neurons_jj);
        uleds_temp.coactZ(kk,jj) = nanmean(temp(:));
        temp = pairwise.postZ(neurons_kk,neurons_jj);
        uleds_temp.postZ(kk,jj) = nanmean(temp(:));
        uleds_temp.post_pre(kk,jj) = uleds_temp.postZ(kk,jj) - uleds_temp.preZ(kk,jj);
    end
end
uLEDcoactivation.uleds_pyramidalCells_lightResponsive = uleds_temp;

% coactivation zscore for uleds and neurons using light responsive
uleds_temp.preZ = [];
uleds_temp.coactZ = [];
uleds_temp.postZ = [];
uleds_temp.post_pre = [];
uleds_temp.rateZmat = [];
for kk = 1:size(lightRespMat,2)
    for jj = 1:size(lightRespMat,2)
        % find neurons
        neurons_kk = find(lightRespMat(:,kk)==1);
        neurons_jj  = find(lightRespMat(:,jj )==1);
        uleds_temp.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

        % coactivations of uleds
        temp = pairwise.preZ(neurons_kk,neurons_jj);
        uleds_temp.preZ(kk,jj) = nanmean(temp(:));
        temp = pairwise.coactZ(neurons_kk,neurons_jj);
        uleds_temp.coactZ(kk,jj) = nanmean(temp(:));
        temp = pairwise.postZ(neurons_kk,neurons_jj);
        uleds_temp.postZ(kk,jj) = nanmean(temp(:));
        uleds_temp.post_pre(kk,jj) = uleds_temp.postZ(kk,jj) - uleds_temp.preZ(kk,jj);
    end
end
uLEDcoactivation.uleds_lightResponsive = uleds_temp;

% coactivation of neurons
temp = zeros(size(lightRespMat,1));
for ii = 1:size(lightRespMat,1)
    for jj = 1:size(lightRespMat,1)
        if any(lightRespMat(ii,:) + lightRespMat(jj,:) == 2) % if both neurons are coactivated at least by one light
            temp(ii,jj) = 1;
        end
    end
end
coactivated_neurons.are_coactivated_pairwise = temp;
coactivated_neurons.are_coactivated_pairs = reshape(temp',[],1);

coactivated_neurons.are_distint_electrode_pairwise = (temp + double(distance.are_distint_electrode_pairwise)) == 2;
coactivated_neurons.are_distint_electrode_pairs = reshape(coactivated_neurons.are_distint_electrode_pairwise',[],1);

coactivated_neurons.are_distint_shank_pairwise = (temp + double(distance.are_distint_shank_pairwise)) == 2;
coactivated_neurons.are_distint_shank_pairs = reshape(coactivated_neurons.are_distint_shank_pairwise',[],1);

coactivated_neurons.are_same_shank_pairwise = (temp + double(distance.are_same_shank_pairwise)) == 2;
coactivated_neurons.are_same_shank_pairs = reshape(coactivated_neurons.are_same_shank_pairwise',[],1);

uLEDcoactivation.coactivated_neurons = coactivated_neurons;


% all neurons
% Coactivation profile of coactivated neurons
uLEDcoactivation.pairwise_coactivated_neurons = mask_unit_as_nan(uLEDcoactivation.pairwise, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_neurons = label_pairs_as_nan(uLEDcoactivation.pairs, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in distinct electrodes
uLEDcoactivation.pairwise_coactivated_neurons_distint_electrode = mask_unit_as_nan(uLEDcoactivation.pairwise_distinct_electrode, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_neurons_distint_electrode = label_pairs_as_nan(uLEDcoactivation.pair_distinct_electrode, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in distinct shank
uLEDcoactivation.pairwise_coactivated_neurons_distint_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_distinct_shank, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_neurons_distint_shank = label_pairs_as_nan(uLEDcoactivation.pair_distinct_shank, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in same shank
uLEDcoactivation.pairwise_coactivated_neurons_same_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_same_shank, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_neurons_same_shank = label_pairs_as_nan(uLEDcoactivation.pair_same_shank, ~coactivated_neurons.are_coactivated_pairs);

% only pyr
% Coactivation profile of coactivated neurons
uLEDcoactivation.pairwise_coactivated_pyramidalCells = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_pyramidalCells = label_pairs_as_nan(uLEDcoactivation.pairs_pyramidalCells, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in distinct electrodes
uLEDcoactivation.pairwise_coactivated_pyramidalCells_distint_electrode = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_distinct_electrode, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_pyramidalCells_distint_electrode = label_pairs_as_nan(uLEDcoactivation.pair_pyramidalCells_distinct_electrode, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in distinct shank
uLEDcoactivation.pairwise_coactivated_pyramidalCells_distint_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_distinct_shank, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_pyramidalCells_distint_shank = label_pairs_as_nan(uLEDcoactivation.pair_pyramidalCells_distinct_shank, ~coactivated_neurons.are_coactivated_pairs);

% Coactivation profile of coactivated neurons in same shank
uLEDcoactivation.pairwise_coactivated_pyramidalCells_same_shank = mask_unit_as_nan(uLEDcoactivation.pairwise_pyramidalCells_same_shank, ~coactivated_neurons.are_coactivated_pairwise);
uLEDcoactivation.pairs_coactivated_pyramidalCells_same_shank = label_pairs_as_nan(uLEDcoactivation.pair_pyramidalCells_same_shank, ~coactivated_neurons.are_coactivated_pairs);

% 1.3 % explained variance
evStats_pairwise = explained_variance(spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~pyramidal_neurons) = [];
evStats_pyramidalCells = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~targetCells) = [];
evStats_pyramidalCells_lightResponsive = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
uLEDcoactivation.evStats_pairwise = evStats_pairwise;
uLEDcoactivation.evStats_pyramidalCells = evStats_pyramidalCells;
uLEDcoactivation.evStats_pyramidalCells_lightResponsive = evStats_pyramidalCells_lightResponsive;

% 1.4. % Groups of neurons
neuronLED = lightRespMat;
neuronLED(isnan(neuronLED) | neuronLED == -1) = 0;
if curate_responses
    neuronLED_curated = curateNeuronLED(neuronLED);
    neuronLED = neuronLED_curated;
end

ledCoactivation = coactivation_matrix;

nNeurons = size(neuronLED,1);
coactivationMatrix_neurons = false(nNeurons);
for i = 1:nNeurons
    for j = i+1:nNeurons
        leds_i = find(neuronLED(i,:));
        leds_j = find(neuronLED(j,:));
        
        activatedTogether = false;
        for li = leds_i
            for lj = leds_j
                if ledCoactivation(li, lj)
                    activatedTogether = true;
                    break
                end
            end
            if activatedTogether, break; end
        end
        
        coactivationMatrix_neurons(i,j) = activatedTogether;
        coactivationMatrix_neurons(j,i) = activatedTogether;  % simetría
    end
end
G = graph(coactivationMatrix_neurons);
artificial_clusters = conncomp(G);
artificial_clusters(degree(G) == 0) = 0;  % marcar como no coactivadas
uLEDcoactivation.clusters.artificial_clusters = artificial_clusters;
uLEDcoactivation.neuronLED_matrix = neuronLED;

% Saving 
if saveMat
    disp('Saving...');
    save([basenameFromBasepath(pwd) '.' save_as '.cellinfo.mat'],'uLEDcoactivation');
end

if doPlot
    
    yaxis_lim = [-10 10];
    % 1 stats
    % all
    figure('Position', [100, 100, 1800, 1000]); % Create and set size in one step
    tiledlayout(3,10)
    nexttile
    groupStats({uLEDcoactivation.pairs.preZ, uLEDcoactivation.pairs.coactZ, uLEDcoactivation.pairs.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylabel('Cofiring (SD)');
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['All neurons'],'FontWeight','normal');
    
    % Pyramidal neurons
    nexttile
    groupStats({uLEDcoactivation.pairs_pyramidalCells.preZ, uLEDcoactivation.pairs_pyramidalCells.coactZ, uLEDcoactivation.pairs_pyramidalCells.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Pyr neurons'],'FontWeight','normal');

    % Light responsive pyr
    nexttile
    groupStats({uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ, uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ, uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Light responsive Pyr'],'FontWeight','normal');

    % Coactive neurons
    nexttile
    groupStats({uLEDcoactivation.pairs_coactivated_neurons.preZ, uLEDcoactivation.pairs_coactivated_neurons.coactZ, uLEDcoactivation.pairs_coactivated_neurons.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Coactive neurons'],'FontWeight','normal');

    % Coact distint electrode
    nexttile
    groupStats({uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.preZ, uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.coactZ, uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Coact different elect'],'FontWeight','normal');

    % Coact distint shank
    nexttile
    groupStats({uLEDcoactivation.pairs_coactivated_neurons_distint_shank.preZ, uLEDcoactivation.pairs_coactivated_neurons_distint_shank.coactZ, uLEDcoactivation.pairs_coactivated_neurons_distint_shank.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Coact different shank'],'FontWeight','normal');

    % Coact same shank
    nexttile
    groupStats({uLEDcoactivation.pairs_coactivated_neurons_same_shank.preZ, uLEDcoactivation.pairs_coactivated_neurons_same_shank.coactZ, uLEDcoactivation.pairs_coactivated_neurons_same_shank.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    ylim(yaxis_lim);
    title(['Coact same shank'],'FontWeight','normal');
    
    yaxis_lim = [-2 2];
    % Abs all
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs.preZ)), log10(abs(uLEDcoactivation.pairs.coactZ)), log10(abs(uLEDcoactivation.pairs.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim(yaxis_lim);
    LogScale('y',10);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    title(['All neurons'],'FontWeight','normal');
    ylabel('abs cofiring (SD)');

    % Abs all
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs_pyramidalCells.preZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim(yaxis_lim);
    LogScale('y',10);
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
    title(['Pyr neurons'],'FontWeight','normal');
    ylabel('abs cofiring (SD)');

    % Abs pyr
    nexttile
    groupStats({log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ)), log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ))},[],'inAxis', true,'plotType','roundPlot','plotData',true, 'repeatedMeasures', true);
    ylim(yaxis_lim);
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

    % correlation coactivated neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons.postZ - uLEDcoactivation.pairs_coactivated_neurons.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons.postZ - uLEDcoactivation.pairs_coactivated_neurons.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % correlation coactivated neurons in different electrodes
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % correlation coactivated neurons in different shanks
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_shank.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_shank.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x, (y),'inAxis', true,'MarkerColor',[.5 .3 .8]);
    xlabel('Coact cofiring (log10 Z)');

    % correlation coactivated neurons in different shanks
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_same_shank.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_same_shank.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_same_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_same_shank.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_same_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_same_shank.preZ) + 1); % Add 1 to avoid log(0)
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

    % 3 positive and negative corr
    % all neurons
    nexttile
    x = sign(uLEDcoactivation.pairs.coactZ) .* log10(abs(uLEDcoactivation.pairs.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) .* log10(abs(uLEDcoactivation.pairs.postZ - uLEDcoactivation.pairs.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');
    
    % pyr neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells.postZ - uLEDcoactivation.pairs_pyramidalCells.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');
    
    % light responsive pyr
    nexttile
    x = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) .* log10(abs(uLEDcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uLEDcoactivation.pairs_pyramidalCells_lightResponsive.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    % coactive neurons
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons.postZ - uLEDcoactivation.pairs_coactivated_neurons.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons.postZ - uLEDcoactivation.pairs_coactivated_neurons.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    % coactive neurons different electrode
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_electrode.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    % coactive neurons different shank
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_shank.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_distint_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_distint_shank.preZ) + 1); % Add 1 to avoid log(0)
    groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.7 .5 1], 'removeOutliers',false);
    groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.3 .1 .6], 'removeOutliers',false, 'labelOffset', 2);
    xlabel('Coact cofiring (log10 Z)');
    ylabel('Post-pre (log10 Z)');

    % coactive neurons same shank
    nexttile
    x = sign(uLEDcoactivation.pairs_coactivated_neurons_same_shank.coactZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_same_shank.coactZ) + 1); % Add 1 to avoid log(0)
    y = sign(uLEDcoactivation.pairs_coactivated_neurons_same_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_same_shank.preZ) .* log10(abs(uLEDcoactivation.pairs_coactivated_neurons_same_shank.postZ - uLEDcoactivation.pairs_coactivated_neurons_same_shank.preZ) + 1); % Add 1 to avoid log(0)
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
    xlabel({['Explain variance (EV vs REV)']; ['all: ' num2str(round(uLEDcoactivation.evStats_pairwise.ev,1)) ' vs ' num2str(round(uLEDcoactivation.evStats_pairwise.rev,1))]; ...
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

function neuronLED_curado = curateNeuronLED(neuronLED)

    [nNeurons, nLEDs] = size(neuronLED);

    data = num2cell(logical(neuronLED));

    f = figure('Name', 'Curado manual de neuron-LED', ...
               'Position', [100 100 100 + 54*nLEDs 100 + 20*nNeurons]);

    % Create table
    t = uitable(f, ...
        'Data', data, ...
        'ColumnEditable', true(1, nLEDs), ...
        'ColumnName', strcat("LED_", string(1:nLEDs)), ...
        'RowName', strcat("N", string(1:nNeurons)), ...
        'Position', [20 60 54*nLEDs 20*nNeurons]);

    % Save button
    uicontrol(f, 'Style', 'pushbutton', ...
        'String', 'Save', ...
        'Position', [20 20 100 30], ...
        'Callback', @(src, event) uiresume(f));

    uiwait(f);

    data_out = get(t, 'Data');
    neuronLED_curado = false(nNeurons, nLEDs);
    for i = 1:nNeurons
        for j = 1:nLEDs
            neuronLED_curado(i,j) = data_out{i,j};
        end
    end

    close(f);
end

