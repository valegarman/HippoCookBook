function plasticity_results = runPlasticityRule_uLED_fcDelta(fc_pre, fc_post, neuronLED_matrix, LED_pulses, varargin)
% plasticity_results = runPlasticityRule_uLED_fcDelta(fc_pre, fc_post, neuronLED_matrix, LED_pulses, varargin)
% 
% Computes changes in functional connectivity (fc_post - fc_pre) for neuron
% pairs coactivated by each LED, based on real pulse times.
% Performs z-scoring against shuffled neuronLED matrices and compares to non-coactivated pairs.
% 
% INPUTS:
%   - fc_pre, fc_post: structs with field .values (NxN) representing connectivity matrices
%   - neuronLED_matrix: NxL logical matrix; neuron i is activated by LED j
%   - LED_pulses: struct with fields:
%         - .timestamps: [M x 2] array with [start, end] of each LED pulse
%         - .code: [M x 1] vector indicating LED index (1 to L) for each pulse
% 
% OPTIONS:
%   - 'min_pairs_per_LED': minimum number of coactivated pairs to include an LED
%   - 'do_plot': whether to plot histogram of delta FCs
%   - 'nShuffles': number of shuffles for null distribution
%   - 'z_thresh': threshold for significant z-score
% 
% OUTPUT:
%   - plasticity_results: struct with fields:
%         - .delta_fc_pairs: vector of delta FC values for coactivated pairs
%         - .pair_idx: [N x 2] indices of neuron pairs
%         - .by_LED: cell array with delta FC values grouped by LED
%         - .null_dist: shuffled delta FC values
%         - .z_scores: z-scored delta FC values
%         - .z_scores_per_pair: per-pair z-scores using shuffled nulls
%         - .coact_mask: logical matrix of coactivated pairs
%         - .non_coact_delta_fc: delta FC for non-coactivated pairs
%         - .significant_pairs: indices of significantly modulated pairs
%         - .significant_z: z-score values for significant pairs
%         - .significant_delta_fc: delta FC values for significant pairs

% Parse inputs
p = inputParser;
addParameter(p, 'min_pairs_per_LED', 3, @isnumeric);
addParameter(p, 'do_plot', true, @islogical);
addParameter(p, 'nShuffles', 1000, @isnumeric);
addParameter(p, 'z_thresh', 2, @isnumeric);
parse(p, varargin{:});

min_pairs_per_LED = p.Results.min_pairs_per_LED;
do_plot = p.Results.do_plot;
nShuffles = p.Results.nShuffles;
z_thresh = p.Results.z_thresh;

N = size(fc_pre.values, 1);
L = size(neuronLED_matrix, 2);
delta_fc = fc_post.values - fc_pre.values;
delta_fc(eye(N) == 1) = NaN;  % remove diagonal

LED_coactivation = false(L, N, N);
for led = 1:L
    neurons_on = find(neuronLED_matrix(:, led));
    for i = 1:length(neurons_on)
        for j = i+1:length(neurons_on)
            LED_coactivation(led, neurons_on(i), neurons_on(j)) = true;
            LED_coactivation(led, neurons_on(j), neurons_on(i)) = true;
        end
    end
end

pulse_delta_fc = [];
pulse_pair_idx = [];
pairs_by_LED = cell(L,1);
coact_mask = false(N,N);

for led = 1:L
    pulses_idx = find(LED_pulses.code == led);
    if isempty(pulses_idx), continue; end

    coact_pairs = squeeze(LED_coactivation(led, :, :));
    [i_idx, j_idx] = find(triu(coact_pairs, 1));
    if length(i_idx) < min_pairs_per_LED, continue; end

    for k = 1:length(i_idx)
        val = delta_fc(i_idx(k), j_idx(k));
        if isnan(val) || val == 0
            continue
        end
        pulse_delta_fc(end+1,1) = val;
        pulse_pair_idx(end+1,:) = [i_idx(k), j_idx(k)];
        pairs_by_LED{led}(end+1) = val;
        coact_mask(i_idx(k), j_idx(k)) = true;
    end
end

% Compute delta FC for non-coactivated pairs
non_coact_mask = triu(~coact_mask, 1);
[i_nc, j_nc] = find(non_coact_mask);
non_coact_delta_fc = arrayfun(@(a,b) delta_fc(a,b), i_nc, j_nc);
non_coact_delta_fc(non_coact_delta_fc == 0) = [];

% Null distribution by shuffling neuronLED_matrix (by column)
nPairs = size(pulse_pair_idx,1);
vals_shuff = nan(nShuffles, nPairs);

for s = 1:nShuffles
    shuffled = false(size(neuronLED_matrix));
    for col = 1:L
        on_idx = find(neuronLED_matrix(:,col));
        shuffled(:,col) = false;
        shuffled(randperm(N, numel(on_idx)), col) = true;
    end

    LED_coactivation_shuff = false(N, N);
    for led = 1:L
        neurons_on = find(shuffled(:, led));
        for i = 1:length(neurons_on)
            for j = i+1:length(neurons_on)
                LED_coactivation_shuff(neurons_on(i), neurons_on(j)) = true;
                LED_coactivation_shuff(neurons_on(j), neurons_on(i)) = true;
            end
        end
    end
    for k = 1:nPairs
        i = pulse_pair_idx(k,1);
        j = pulse_pair_idx(k,2);
        if LED_coactivation_shuff(i,j) && delta_fc(i,j) ~= 0
            vals_shuff(s,k) = delta_fc(i,j);
        end
    end
end

% Z-scores per pair (robust check)
z_scores_per_pair = nan(nPairs,1);
for k = 1:nPairs
    null_vals = vals_shuff(:,k);
    null_vals = null_vals(~isnan(null_vals));
    if length(null_vals) >= 10 && std(null_vals) > 0
        z_scores_per_pair(k) = (pulse_delta_fc(k) - mean(null_vals)) / std(null_vals);
    else
        z_scores_per_pair(k) = NaN;
    end
end

% === Diagnóstico visual de distribuciones null ===
figure;
for kk = 1:min(3,nPairs)
    subplot(1,3,kk)
    null_vals = vals_shuff(:,kk);
    null_vals = null_vals(~isnan(null_vals));
    histogram(null_vals, 30);
    hold on;
    xline(pulse_delta_fc(kk), 'r', 'LineWidth', 2);
    title(sprintf('Par %d vs %d', pulse_pair_idx(kk,1), pulse_pair_idx(kk,2)));
    xlabel('\Delta FC null'); ylabel('count');
end

valid_counts = sum(~isnan(vals_shuff), 1);
figure;
histogram(valid_counts);
xlabel('# null values per pair'); ylabel('# pairs');
title('Número de valores válidos en null por par');

figure;
histogram(z_scores_per_pair, 40);
xlabel('Z-score per pair'); ylabel('# pairs');
title('Z-score distribution (per-pair null)');

% Global z-score (legacy)
mu_null = mean(vals_shuff(:), 'omitnan');
sigma_null = std(vals_shuff(:), 'omitnan');
z_scores = (pulse_delta_fc - mu_null) / sigma_null;

% Identify significant pairs
significant_mask = z_scores_per_pair > z_thresh;

% Output
plasticity_results.delta_fc_pairs = pulse_delta_fc;
plasticity_results.pair_idx = pulse_pair_idx;
plasticity_results.by_LED = pairs_by_LED;
plasticity_results.null_dist = vals_shuff;
plasticity_results.z_scores = z_scores;
plasticity_results.z_scores_per_pair = z_scores_per_pair;
plasticity_results.coact_mask = coact_mask;
plasticity_results.non_coact_delta_fc = non_coact_delta_fc;
plasticity_results.significant_pairs = pulse_pair_idx(significant_mask,:);
plasticity_results.significant_z = z_scores_per_pair(significant_mask);
plasticity_results.significant_delta_fc = pulse_delta_fc(significant_mask);

if do_plot
    figure;
    histogram(non_coact_delta_fc, 40, 'FaceColor', [0.85 0.4 0.3], 'FaceAlpha', 0.6, 'DisplayName','non-coactivated'); hold on;
    histogram(pulse_delta_fc, 40, 'FaceColor', [0.2 0.5 0.9], 'FaceAlpha', 0.6, 'DisplayName','coactivated');
    xlabel('\Delta FC (post - pre)'); ylabel('# pairs'); title('Delta FC distributions'); legend;
end
end