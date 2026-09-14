
function [functional_connectivity] = getFunctionalConnectivity(varargin)
% Compute functional connectivity between neuron families as in Valero, Abad et
% al, 2025

% INPUTS

% OUTPUTS

% MV-NCL2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[]);
addParameter(p,'method','te',@ischar); % options are mi, glm, te
addParameter(p,'restrict',[],@isnumeric); % ints to restrict to
addParameter(p,'dt',0.02,@isnumeric); % 
addParameter(p,'k',3,@isnumeric); % history length in bins, for transfer entropy
addParameter(p,'n_shuffles',1000,@isnumeric); 
addParameter(p,'doPlot',true,@islogical); 
addParameter(p,'alpha',0.01,@isnumeric); 
addParameter(p,'save_as','functionalConnectivity',@ischar); 
addParameter(p,'savemat',true,@islogical); 
addParameter(p,'ccg_win',[-10 10], @isnumeric); 
        
parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
method = p.Results.method;
restrict = p.Results.restrict;
dt = p.Results.dt;
k = p.Results.k;
n_shuffles = p.Results.n_shuffles;
doPlot = p.Results.doPlot;
alpha = p.Results.alpha;
save_as = p.Results.save_as;
savemat = p.Results.savemat;
ccg_win = p.Results.ccg_win;

% 
prevPath = pwd;
cd(basepath)

if isempty(spikes)
    spikes = loadSpikes;
end

if ~isempty(restrict)
    for ii = 1:spikes.numcells
        spikes.times{ii} = Restrict(spikes.times{ii}, restrict);
    end
end

disp('Estimating functional connectivity...');
spikemat = bz_SpktToSpkmat(spikes,'dt',dt);
functional_connectivity.values = [];
functional_connectivity.p = [];
functional_connectivity.method = method;
functional_connectivity.n_shuffles = n_shuffles;
switch lower(method)
    case 'mi'
        disp('... mi...');
        % MI
        nCols = size(spikemat.data, 2);
        nTime = size(spikemat.data, 1);
        n_perm = n_shuffles;
        
        mi_mat = nan(nCols, nCols);
        mi_p = nan(nCols, nCols);
        
        spikemat_data = int32(spikemat.data);  % cast once to int32
        
        % Precompute permutation indices (can be reused)
        perm_idx = zeros(nTime, n_perm, 'uint32');
        for p = 1:n_perm
            perm_idx(:, p) = randperm(nTime);
        end
        
        for jj = 1:nCols
            fprintf('Computing MI: neuron %d of %d\n', jj, nCols);
            for kk = jj:nCols  % use upper triangle only
                x = spikemat_data(:, jj);
                y = spikemat_data(:, kk);
        
                % Observed MI
                mi_obs = nmi(x', y');
                mi_mat(jj, kk) = mi_obs;
                mi_mat(kk, jj) = mi_obs;  % symmetry
        
                % Permutation null
                mi_null = nan(n_perm, 1);
                for p = 1:n_perm
                    y_perm = y(perm_idx(:, p));
                    mi_null(p) = nmi(x', y_perm');
                end
        
                % p-value
                p_val = mean(mi_null >= mi_obs);
                mi_p(jj, kk) = p_val;
                mi_p(kk, jj) = p_val;
            end
        end
        
        functional_connectivity.values = mi_mat;
        functional_connectivity.p = mi_p;

    case 'glm'
        % GLM
        disp('... glm...');
        nCols = size(spikemat.data, 2);
        T = size(spikemat.data, 1);
        warning('off', 'stats:glmfit:IterationLimit');
        warning('off', 'stats:glmfit:IllConditioned');
        
        % Step 1: Select 10 random target neurons
        rng(1);  % for reproducibility
        subset_neurons = randsample(nCols, 10);
        lambda_values = nan(length(subset_neurons), 1);
        
        for i = 1:length(subset_neurons)
            jj = subset_neurons(i);
            predictorCells = ones(nCols, 1);
            predictorCells(jj) = 0;
        
            X = spikemat.data(:, find(predictorCells));
            Y = spikemat.data(:, find(predictorCells == 0));
        
            try
                [~, FitInfo] = lassoglm(X, Y, 'poisson', 'CV', 5);
                idx = FitInfo.IndexMinDeviance;
                lambda_values(i) = FitInfo.Lambda(idx);
            catch ME
                warning('CV failed for neuron %d: %s', jj, ME.message);
            end
        end
        
        % Step 2: Use the median Lambda for all neurons
        Lambda_best = median(lambda_values(~isnan(lambda_values)));
        fprintf('Using Lambda = %.6f for all neurons\n', Lambda_best);
        
        glm_mat = nan(nCols, nCols);
        glm_p = nan(nCols, nCols);
        
        fprintf('\n');
        for jj = 1:nCols
            fprintf('\rComputing GLM: neuron %d of %d', jj, nCols);
            predictorCells = ones(nCols, 1);
            predictorCells(jj) = 0;
            opts = statset('glmfit');
            opts.MaxIter = 1000;

            X = spikemat.data(:, find(predictorCells));
            Y = spikemat.data(:, find(predictorCells == 0));
        
            try
                [B, FitInfo] = lassoglm(X, Y, 'poisson', 'Lambda', Lambda_best);
                selected = find(B ~= 0);
        
                if ~isempty(selected)
                    X_sel = zscore(X(:, selected));
                    keep = std(X_sel) > 1e-8;
                    X_sel = X_sel(:, keep);
                   
                    % Fit null model (intercept only)
                    [~, dev_null] = glmfit(ones(size(Y)), Y, 'poisson', 'Options', opts);
                    
                    % Fit full model
                    [b_glm, dev, stats] = glmfit(X_sel, Y, 'poisson', 'Options', opts);

                    % Compute pseudo R²
                    glm_pseudoR2(jj) = 1 - (dev / dev_null);
        
                    selectedIndices = find(predictorCells);
                    selected_orig = selectedIndices(selected);
        
                    glm_mat(selected_orig, jj) = b_glm(2:end);
                    glm_p(selected_orig, jj) = stats.p(2:end);

                    % Store model deviance for target neuron jj
                    glm_deviance(jj) = dev;
                end
            
            catch ME
                warning('GLM failed for neuron %d: %s', jj, ME.message);
            end
        end
        
        functional_connectivity.values = glm_mat;
        functional_connectivity.p = glm_p;
        functional_connectivity.deviance = glm_deviance;
        functional_connectivity.pseudoR2  = glm_pseudoR2;
        warning('on', 'stats:glmfit:IterationLimit');
        warning('on', 'stats:glmfit:IllConditioned');
    
    case 'glm_conditional'
    disp('... glm with conditional population covariate...');
    nCols = size(spikemat.data, 2);
    glm_cond_mat = nan(nCols, nCols);
    glm_cond_p = nan(nCols, nCols);
    glm_cond_pseudoR2 = nan(nCols, nCols);
    glm_cond_deviance = nan(nCols, nCols);
    warning('off', 'stats:glmfit:IterationLimit');
    warning('off', 'stats:glmfit:IllConditioned');

    % Fit with higher iteration limit
    opts = statset('glmfit');
    opts.MaxIter = 1000;
    fprintf('\n');
    for j = 1:nCols
        Y = spikemat.data(:, j);  % target neuron
        fprintf('\rComputing GLM: neuron %d of %d', j, nCols);
        for i = 1:nCols
            if i == j
                continue;  % skip self-pairs
            end

            % Predictor 1: activity of source neuron i
            x = spikemat.data(:, i);

            % Predictor 2: summed population activity excluding i and j
            exclude = true(1, nCols);
            exclude([i j]) = false;
            z = sum(spikemat.data(:, exclude), 2);

            % Design matrix: predictors [x, z]
            X = zscore([x z]);  % Apply zscore here
            keep = std(X) > 1e-8;
            X = X(:, keep);

            try
                % Fit null model (intercept only)
                try
                    [~, dev_null] = glmfit(ones(size(Y)), Y, 'poisson', 'Options', opts);
                catch
                    dev_null = NaN;
                    warning('Null model failed for neuron %d', jj);
                end

                % Fit full model with predictors x and z
                [b, dev, stats] = glmfit(X, Y, 'poisson', 'Options', opts);

                % Store results
                glm_cond_mat(i, j) = b(2);              % β: effect of i on j
                glm_cond_p(i, j) = stats.p(2);          % p-value for β
                glm_cond_pseudoR2(i, j) = 1 - (dev / dev_null);
                glm_cond_deviance(i, j) = dev;

            catch ME
                fprintf('\n');
                warning('Conditional GLM failed for i = %d → j = %d: %s', i, j, ME.message);
                fprintf('\n');
            end
        end
    end
    fprintf('\n');
    % Store in output struct
    functional_connectivity.values = glm_cond_mat;
    functional_connectivity.p = glm_cond_p;
    functional_connectivity.pseudoR2 = glm_cond_pseudoR2;
    functional_connectivity.deviance = glm_cond_deviance;
    functional_connectivity.direction = 'i_to_j';
    warning('on', 'stats:glmfit:IterationLimit');
    warning('on', 'stats:glmfit:IllConditioned');
    
    case 'te'

        % Transfer entropy
        disp('... transfer entropy...');
        nCols = size(spikemat.data, 2);
        te_mat =  nan(size(spikemat.data, 2), size(spikemat.data, 2));
        te_p =  nan(size(spikemat.data, 2), size(spikemat.data, 2));
        for jj = 1:nCols
            fprintf('Computing TE: neuron %d of %d\n', jj, nCols);
            source = double(spikemat.data(:, jj));
            for kk =  1:nCols
                if jj == kk
                    continue
                end
                target = double(spikemat.data(:, kk));
                temp = compute_TE_counts(source, target, k, n_shuffles);
                % Use raw TE or corrected TE depending on what you want
                te_mat(jj, kk) = temp.value;
                % te_mat(jj, kk) = temp.corrected;

                te_p(jj, kk) = temp.p;
            end
        end

        functional_connectivity.values = te_mat;
        functional_connectivity.p = te_p;

    case 'ccg'
       disp('... CCG method using compiled CCGHeart...');

        % Parameters
        binSize = 0.001;       % in seconds (e.g., 1 ms bins)
        winSize = 0.1;         % total window size (±50 ms)
        
        % Setup
        nCols = spikes.numcells;
        session = loadSession;
        Fs = 1 / session.extracellular.sr;
        
        % Compute all CCGs
        [allCcg, t_ccg] = CCG(spikes.times, [], ...
                              'binSize', binSize, ...
                              'duration', winSize, ...
                              'Fs', Fs, ...
                              'normtype', 'counts');
        
        % Get center bin
        center_bin = ceil((winSize / binSize + 1) / 2);
        
        % Apply custom window
        offset_bins = ccg_win(1):ccg_win(2);
        selected_bins = center_bin + offset_bins;
        selected_bins = selected_bins(selected_bins > 0 & selected_bins <= size(allCcg,1));
        
        % Average across selected bins (raw coactivation)
        ccg_raw = squeeze(mean(allCcg(selected_bins, :, :), 1));  % NxN
        
        % Normalize by expected coactivation
        non_empty_times = spikes.times(~cellfun(@isempty, spikes.times));
        max_time = max(cellfun(@max, non_empty_times));
        min_time = min(cellfun(@min, non_empty_times));
        total_time = max_time - min_time;
        fr = cellfun(@(x) numel(x) / total_time, spikes.times);  % Hz
        expected_coactivation = fr(:) * fr(:)' * (binSize * length(selected_bins));
        ccg_norm = ccg_raw ./ expected_coactivation;
        
        % Shuffling (circular)
        n_shuffles = p.Results.n_shuffles;
        ccg_null = zeros(n_shuffles, nCols, nCols);
        rng(1);  % reproducibility
        
        for s = 1:n_shuffles
            shuffled_times = cellfun(@(x) mod(x + rand * total_time, total_time), ...
                                     spikes.times, 'UniformOutput', false);
            [ccg_s, ~] = CCG(shuffled_times, [], ...
                             'binSize', binSize, ...
                             'duration', winSize, ...
                             'Fs', Fs, ...
                             'normtype', 'counts');
            ccg_null(s,:,:) = squeeze(mean(ccg_s(selected_bins, :, :), 1));
        end
        
        % Null distribution statistics
        mu_null = squeeze(mean(ccg_null, 1));
        std_null = squeeze(std(ccg_null, 0, 1));
        z_mat = (ccg_raw - mu_null) ./ std_null;
        p_mat = squeeze(mean(ccg_null >= reshape(ccg_raw, [1 nCols nCols]), 1));  % one-tailed
        DI_z = (z_mat - z_mat') ./ (z_mat + z_mat');
        
        % Assign outputs
        functional_connectivity.values = ccg_norm;
        functional_connectivity.raw = ccg_raw;
        functional_connectivity.z = z_mat;
        functional_connectivity.p = p_mat;
        functional_connectivity.ccg_window = ccg_win;
        functional_connectivity.directionality_indexZ = DI_z;
        functional_connectivity.direction = '';
    otherwise
        error('Method do not recognized...');
end

%% Classify cell types
[cell_types, cell_classification_stats, cell_metrics, cell_subtypes] = cellTypeClassifier('modelType','hippocampus5');
% Unique classes (excluding 'Noisy_unit')
cell_classes = unique(cell_types);
cell_classes(strcmp(cell_classes, 'Noisy_unit')) = [];  % remove noisy
cell_classes(strcmp(cell_classes, 'Undetermined')) = [];  % remove undetermined


n_classes = length(cell_classes);
class_conn_avg = nan(n_classes, n_classes);  % rows: target, cols: source
class_conn_p = nan(n_classes, n_classes);  % average p-value per class pair

% 
funcContSig = functional_connectivity.values;
funcContSig(functional_connectivity.p > alpha) = NaN;  % 
for tgt = 1:n_classes
    for src = 1:n_classes
        % Logical index of cells in each class
        tgt_idx = strcmp(cell_types, cell_classes{tgt});
        src_idx = strcmp(cell_types, cell_classes{src});

        % Extract relevant submatrix
        conn_block = funcContSig(tgt_idx, src_idx);
        
        % Exclude NaNs and zeros if desired
        conn_vals = conn_block(~isnan(conn_block));

        % Extract the block of p-values
        p_block = functional_connectivity.p(tgt_idx, src_idx);
        p_vals = p_block(~isnan(p_block));  % remove NaNs

        if ~isempty(conn_vals)
            class_conn_avg(tgt, src) = nanmean(conn_vals);
        end

        if ~isempty(p_vals)
            class_conn_p(tgt, src) = nanmean(p_vals);
        end
    end
end

functional_connectivity.cell_families_conn_avg = class_conn_avg;
functional_connectivity.cell_families_conn_p = class_conn_p;
functional_connectivity.cell_families = cell_classes;

if savemat
    disp('Saving results...');
    save([basenameFromBasepath(pwd) '.' save_as '.cellinfo.mat'], 'functional_connectivity');
end

%% Plots
if doPlot
    switch lower(method)
        case 'mi'
            figure
            subplot(1,2,1)
            plotCorrMat(functional_connectivity.cell_families_conn_avg, functional_connectivity.cell_families_conn_p,'minPvalue', 1E-20,'maxPvalue',0.1, 'area_factor', 7,...
            'Y_variablesNames', functional_connectivity.cell_families,...
            'X_variablesNames', functional_connectivity.cell_families, 'p_value_threshold', 0.0001, 'inAxis', true);   
            caxis([0 .02]);
            colormap(flip(gray));
            axis square
            set(gca, 'TickLabelInterpreter', 'none');
            
            x_axis = functional_connectivity.cell_families;
            avg = functional_connectivity.cell_families_conn_avg;
            
            % Optional: remove self-connections (comment this out if you want to include them)
            % avg(1:size(avg,1)+1:end) = 0;
            
            % Symmetrize matrix
            symMatrix = (avg + avg') / 2;
            colors = getColors;

            % Build and plot graph
            subplot(1,2,2)
            G = graph(symMatrix);
            h = plot(G, 'Layout', 'layered');  % 'circle' layout looks cleaner
            
            % Use sqrt scaling for better visibility
            weights = G.Edges.Weight;
            scaled_weights = sqrt(weights);
            scaled_weights = scaled_weights / max(scaled_weights);  % Normalize
            
            % Line width based on scaled weights
            h.LineWidth = 3; % 2 + 8 * scaled_weights;  % Range: [2,10]
            
            % Edge color from grayscale colormap
            num_colors = 256;
            cmap = flip(gray(num_colors));
            indices = round(scaled_weights * (num_colors - 1)) + 1;
            h.EdgeColor = cmap(indices, :);
            
            % Match node labels
            h.NodeLabel = x_axis;
            
            % Match node colors — adjust this to fit actual order of classes
            % For example, if x_axis = {'CAMK2_DEEP', 'CAMK2_SUP', 'PV+', 'SST+'};
            h.NodeColor = [colors.deep; colors.sup; colors.pv; colors.sst];
            
            % Aesthetics
            h.MarkerSize = 8;
            axis off

            mkdir('SummaryFigures');
            exportgraphics(gcf,['SummaryFigures\' save_as '_' method '.png']);
        
        case {'glm', 'glm_conditional', 'te'}
            figure
            subplot(1,2,1)
            plotCorrMat(functional_connectivity.cell_families_conn_avg, functional_connectivity.cell_families_conn_p,'minPvalue', 1E-20,'maxPvalue',0.1, 'area_factor', 7,...
            'Y_variablesNames', functional_connectivity.cell_families,...
            'X_variablesNames', functional_connectivity.cell_families, 'p_value_threshold', 0.0001, 'inAxis', true);   
            % caxis([-0.05 .05]);
            axis square
            set(gca, 'TickLabelInterpreter', 'none');
            
            x_axis = functional_connectivity.cell_families;
            avg = functional_connectivity.cell_families_conn_avg;
            avg(isnan(avg)) = 0;
            colors = getColors(functional_connectivity.cell_families);
            G = digraph(avg);
            subplot(1,2,2)
            h = plot(G, 'Layout', 'layered'); % force circle
            weights = G.Edges.Weight;
            maxAbsWeight = max(abs(weights));
            normalizedWeights = (weights + maxAbsWeight) / (2 * maxAbsWeight); % Centered around 0
            normalizedWeights(normalizedWeights<0) = 0; normalizedWeights(normalizedWeights>1) = 1;
            cmap = colormap(flip(brewermap(100,'RdYlBu'))); % 
            edgeColors = cmap(round(normalizedWeights * 99) + 1, :);
            h.EdgeColor = edgeColors;
            h.LineWidth = 3; % 5 * abs(weights) / max(abs(weights));
            h.NodeColor = colors; % ''; %  x_axis'
            % h.NodeLabel = functional_connectivity.cell_families; % ''; %  x_axis'
            h.NodeLabel = '';
            h.EdgeAlpha = 0.8;
            h.MarkerSize = 8;
            h.ArrowSize = 15;
            axis off

            mkdir('SummaryFigures');
            exportgraphics(gcf,['SummaryFigures\' save_as '_' method '.png']);
        case 'ccg'
            figure
            subplot(1,2,1)
            plotCorrMat(functional_connectivity.cell_families_conn_avg, functional_connectivity.cell_families_conn_p,...
                'minPvalue', 1E-20, 'maxPvalue', 0.1, 'area_factor', 7,...
                'Y_variablesNames', functional_connectivity.cell_families,...
                'X_variablesNames', functional_connectivity.cell_families,...
                'p_value_threshold', 0.01, 'inAxis', true);   
            caxis([0 5]);  % ajusta según tu rango
            axis square
            set(gca, 'TickLabelInterpreter', 'none');
            
            % Red como grafo
            x_axis = functional_connectivity.cell_families;
            avg = functional_connectivity.cell_families_conn_avg;
            z_scores = functional_connectivity.cell_families_conn_avg;  % 

            colors = getColors(x_axis);
            G = digraph(avg);  %

            subplot(1,2,2)
            h = plot(G, 'Layout', 'layered');

            % Grosor según fuerza de conexión
            weights = G.Edges.Weight;
            weights(~isfinite(weights)) = 0;
            maxW = max(weights);
            if maxW > 0
                scaled_weights = sqrt(weights) / sqrt(maxW);
            else
                scaled_weights = zeros(size(weights));  % todo 0 si no hay conexiones
            end
            h.LineWidth = 1 + 4 * scaled_weights;

            % Color según z-score promedio por clase (puedes personalizar)
            num_colors = 256;
            cmap = parula(num_colors);
            normalized_z = (weights - min(weights)) / (max(weights) - min(weights));
            indices = round(normalized_z * (num_colors - 1)) + 1;
            h.EdgeColor = cmap(indices, :);

            h.NodeColor = colors;
            h.NodeLabel = x_axis;
            h.MarkerSize = 8;
            axis off

            mkdir('SummaryFigures');
            exportgraphics(gcf,['SummaryFigures\' save_as '_' method '.png']);
        otherwise
            warning('Option not implemented yet...');
    end
           
end


cd(prevPath);
end

function TE = compute_TE_counts(source, target, k, n_shuffles)
% COMPUTE_TE_COUNTS Compute transfer entropy from source to target.
%
% TE = I(target_future ; source_past | target_past)
%
% Inputs:
%   source, target : spike count vectors, same length
%   k              : number of past bins
%   n_shuffles     : number of source shuffles
%
% Outputs:
%   TE.value       : raw transfer entropy, in bits
%   TE.null        : shuffled null distribution
%   TE.p           : one-sided shuffle p-value
%   TE.corrected   : TE.value - mean(TE.null)

    if nargin < 4 || isempty(n_shuffles)
        n_shuffles = 100;
    end

    source = source(:);
    target = target(:);

    assert(numel(source) == numel(target), ...
        'source and target must have the same length');

    N = numel(source);
    assert(N > k + 1, 'Time series is too short for this k');

    % Build aligned past/future states manually
    nSamples = N - k;

    X_past = zeros(nSamples, k);
    Y_past = zeros(nSamples, k);
    Y_future = zeros(nSamples, 1);

    for t = 1:nSamples
        idxPast = t:(t + k - 1);

        X_past(t, :) = source(idxPast);
        Y_past(t, :) = target(idxPast);
        Y_future(t)  = target(t + k);
    end

    % Raw TE: I(target_future ; source_past | target_past)
    TE.value = conditional_mutual_info(Y_future, X_past, Y_past);

    % Null distribution by shuffling source
    TE.null = zeros(n_shuffles, 1);

    for s = 1:n_shuffles
        shuffled_source = source(randperm(N));

        X_past_shuff = zeros(nSamples, k);

        for t = 1:nSamples
            idxPast = t:(t + k - 1);
            X_past_shuff(t, :) = shuffled_source(idxPast);
        end

        TE.null(s) = conditional_mutual_info(Y_future, X_past_shuff, Y_past);
    end

    % One-sided p-value
    TE.p = mean(TE.null >= TE.value);

    % Bias-corrected / shuffle-corrected TE
    TE.corrected = TE.value - mean(TE.null);
end

function I = conditional_mutual_info(Z, Y, X)
% CONDITIONAL_MUTUAL_INFO Compute I(Z ; Y | X), in bits.
%
% Inputs:
%   Z : vector, variable of interest
%   Y : vector or matrix, informative variable
%   X : vector or matrix, conditioning variable
%
% Output:
%   I : conditional mutual information in bits
%
% This computes:
%
%   I(Z;Y|X) = sum p(z,y,x) log2( p(z,y,x) p(x) / (p(z,x) p(y,x)) )

    Z = Z(:);

    if isvector(Y)
        Y = Y(:);
    end

    if isvector(X)
        X = X(:);
    end

    n = numel(Z);

    assert(size(Y, 1) == n, 'Z and Y must have the same number of samples');
    assert(size(X, 1) == n, 'Z and X must have the same number of samples');

    % Convert states to unique integer IDs
    z_id = state_id(Z);
    y_id = state_id(Y);
    x_id = state_id(X);

    xyz_states = [z_id, y_id, x_id];
    zx_states  = [z_id, x_id];
    yx_states  = [y_id, x_id];

    [~, ~, xyz_id] = unique(xyz_states, 'rows');
    [~, ~, zx_id]  = unique(zx_states,  'rows');
    [~, ~, yx_id]  = unique(yx_states,  'rows');
    [~, ~, x_id2]  = unique(x_id,       'rows');

    % Empirical probabilities
    p_xyz = accumarray(xyz_id, 1) / n;
    p_zx  = accumarray(zx_id,  1) / n;
    p_yx  = accumarray(yx_id,  1) / n;
    p_x   = accumarray(x_id2,  1) / n;

    % Sum over unique xyz states
    I = 0;

    for state = 1:numel(p_xyz)

        idx = find(xyz_id == state, 1, 'first');

        p1 = p_xyz(state);
        p2 = p_zx(zx_id(idx));
        p3 = p_yx(yx_id(idx));
        p4 = p_x(x_id2(idx));

        if p1 > 0 && p2 > 0 && p3 > 0 && p4 > 0
            I = I + p1 * log2((p1 * p4) / (p2 * p3));
        end
    end

    % Remove tiny numerical negatives
    if I < 0 && abs(I) < 1e-12
        I = 0;
    end
end

function id = state_id(A)
% STATE_ID Convert vector or matrix rows into integer state IDs.

    if isvector(A)
        A = A(:);
    end

    [~, ~, id] = unique(A, 'rows');
end

function ids = vector_to_id(vec)
    ids = vec(:) + 1;
end

function ids = matrix_to_id(mat)
    % Each row is a history vector (counts or binary)
    % Convert each row to unique ID using base conversion
    max_vals = max(mat, [], 1);
    base_vec = cumprod([1, max_vals(1:end-1)+1]);  % dynamic base per column
    ids = mat * base_vec(:) + 1;
end

function h = joint_hist3(a, b)
    max_a = max(a);
    max_b = max(b);
    h = accumarray([a(:), b(:)], 1, [max_a, max_b]);
    h = h / sum(h(:));  % normalize
end