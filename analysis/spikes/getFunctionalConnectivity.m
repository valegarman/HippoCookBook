
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
addParameter(p,'dt',0.1,@isnumeric); % 
addParameter(p,'k',3,@isnumeric); % history length in bins, for transfer entropy
addParameter(p,'n_shuffles',100,@isnumeric); 
addParameter(p,'doPlot',true,@islogical); 
addParameter(p,'alpha',0.01,@isnumeric); 
addParameter(p,'save_as','functionalConnectivity',@ischar); 
addParameter(p,'savemat',true,@islogical); 
        
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
        
        for jj = 1:nCols
            fprintf('Computing GLM: neuron %d of %d\n', jj, nCols);
            predictorCells = ones(nCols, 1);
            predictorCells(jj) = 0;
        
            X = spikemat.data(:, find(predictorCells));
            Y = spikemat.data(:, find(predictorCells == 0));
        
            try
                [B, FitInfo] = lassoglm(X, Y, 'poisson', 'Lambda', Lambda_best);
                selected = find(B ~= 0);
        
                if ~isempty(selected)
                    X_sel = X(:, selected);
                    [b_glm, dev, stats] = glmfit(X_sel, Y, 'poisson');
        
                    selectedIndices = find(predictorCells);
                    selected_orig = selectedIndices(selected);
        
                    glm_mat(jj, selected_orig) = b_glm(2:end);
                    glm_p(jj, selected_orig) = stats.p(2:end);
                end
        
            catch ME
                warning('GLM failed for neuron %d: %s', jj, ME.message);
            end
        end

        functional_connectivity.values = glm_mat;
        functional_connectivity.p = glm_p;
    case 'te'

        % Transfer entropy
        disp('... transfer entropy...');
        nCols = size(spikemat.data, 2);
        te_mat =  nan(size(spikemat.data, 2), size(spikemat.data, 2));
        te_p =  nan(size(spikemat.data, 2), size(spikemat.data, 2));
        for jj = 1:nCols
            fprintf('Computing TE: neuron %d of %d\n', jj, nCols);
            for kk =  1:nCols
                temp = compute_TE_counts(double(spikemat.data(1:end, jj)), double(spikemat.data(1:end, kk)), k, n_shuffles);  % k=2, 100 shuffles
                te_mat(jj,kk) = temp.value;
                te_p(jj,kk) = temp.p;
            end
        end

        functional_connectivity.values = te_mat;
        functional_connectivity.p = te_p;

    case 'ccg'
        disp('... CCG method using compiled CCGHeart...');
        
        % Preparación
        nCols = spikes.numcells;
        binSize = 0.001;     % en segundos
        winSize = 0.1;       % ventana total (±50 ms)
        center_bin = ceil((winSize / binSize + 1) / 2);
        session = loadSession;
        Fs = 1 / session.extracellular.sr;

        % Calcula todos los CCGs simultáneamente
        [allCcg, t_ccg] = CCG(spikes.times, [], ...
                              'binSize', binSize, ...
                              'duration', winSize, ...
                              'Fs', Fs, ...
                              'normtype', 'counts');

        % Extrae el valor del bin central (lag ≈ 0)
        ccg_raw = squeeze(allCcg(center_bin, :, :));  % NxN

        % Normalización por tasa de disparo esperada
        non_empty_times = spikes.times(~cellfun(@isempty, spikes.times));
        max_time = max(cellfun(@max, non_empty_times));
        min_time = min(cellfun(@min, non_empty_times));
        total_time = max_time - min_time;
        fr = cellfun(@(x) numel(x) / total_time, spikes.times);  % Hz
        expected_coactivation = fr(:) * fr(:)' * winSize;
        ccg_norm = ccg_raw ./ expected_coactivation;

        % Shuffling (desplazamiento circular)
        n_shuffles = p.Results.n_shuffles;
        ccg_null = zeros(n_shuffles, nCols, nCols);
        rng(1);  % para reproducibilidad

        for s = 1:n_shuffles
            shuffled_times = cellfun(@(x) mod(x + rand * total_time, total_time), ...
                                     spikes.times, 'UniformOutput', false);
            [ccg_s, ~] = CCG(shuffled_times, [], ...
                             'binSize', binSize, ...
                             'duration', winSize, ...
                             'Fs', Fs, ...
                             'normtype', 'counts');
            ccg_null(s,:,:) = squeeze(ccg_s(center_bin, :, :));
        end

        mu_null = squeeze(mean(ccg_null, 1));
        std_null = squeeze(std(ccg_null, 0, 1));
        z_mat = (ccg_raw - mu_null) ./ std_null;
        p_mat = squeeze(mean(ccg_null >= reshape(ccg_raw, [1 nCols nCols]), 1));  % one-tailed

        % Asignar resultados
        functional_connectivity.values = ccg_norm;
        functional_connectivity.raw = ccg_raw;
        functional_connectivity.z = z_mat;
        functional_connectivity.p = p_mat;
    otherwise
        error('Method do not recognized...');
end

%% Classify cell types
[cell_types, cell_classification_stats, cell_metrics, cell_subtypes] = cellTypeClassifier('modelType','hippocampus5');
% Unique classes (excluding 'Noisy_unit')
cell_classes = unique(cell_subtypes);
cell_classes(strcmp(cell_classes, 'Noisy_unit')) = [];  % remove noisy

n_classes = length(cell_classes);
class_conn_avg = nan(n_classes, n_classes);  % rows: target, cols: source
class_conn_p = nan(n_classes, n_classes);  % average p-value per class pair

% 
funcContSig = functional_connectivity.values;
funcContSig(functional_connectivity.p > alpha) = NaN;  % 
for tgt = 1:n_classes
    for src = 1:n_classes
        % Logical index of cells in each class
        tgt_idx = strcmp(cell_subtypes, cell_classes{tgt});
        src_idx = strcmp(cell_subtypes, cell_classes{src});

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
            class_conn_p(tgt, src) = mean(p_vals);
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
        case 'glm'
            figure
            subplot(1,2,1)
            plotCorrMat(functional_connectivity.cell_families_conn_avg, functional_connectivity.cell_families_conn_p,'minPvalue', 1E-20,'maxPvalue',0.1, 'area_factor', 7,...
            'Y_variablesNames', functional_connectivity.cell_families,...
            'X_variablesNames', functional_connectivity.cell_families, 'p_value_threshold', 0.0001, 'inAxis', true);   
            caxis([-0.05 .05]);
            axis square
            set(gca, 'TickLabelInterpreter', 'none');
            
            x_axis = functional_connectivity.cell_families;
            avg = functional_connectivity.cell_families_conn_avg;
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
            h.NodeLabel = functional_connectivity.cell_families; % ''; %  x_axis'
            h.EdgeAlpha = 0.8;
            h.MarkerSize = 8;
            h.ArrowSize = 15;
            axis off

            mkdir('SummaryFigures');
            exportgraphics(gcf,['SummaryFigures\' save_as '_' method '.png']);
        case 'te'
            figure
            subplot(1,2,1)
            plotCorrMat(functional_connectivity.cell_families_conn_avg, functional_connectivity.cell_families_conn_p,'minPvalue', 1E-20,'maxPvalue',0.1, 'area_factor', 20,...
            'Y_variablesNames', functional_connectivity.cell_families,...
            'X_variablesNames', functional_connectivity.cell_families, 'p_value_threshold', 0.0001, 'inAxis', true);   
            caxis([-0.05 .05]);
            axis square
            set(gca, 'TickLabelInterpreter', 'none');
            
            x_axis = functional_connectivity.cell_families;
            avg = functional_connectivity.cell_families_conn_avg;
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
            h.NodeLabel = functional_connectivity.cell_families; % ''; %  x_axis'
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
            z_scores = functional_connectivity.cell_families_conn_avg;  % puedes poner otra cosa si lo agrupas

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
% Compute TE from 'source' to 'target' with spike count vectors
% Inputs:
%   - source, target: spike count vectors (same length)
%   - k: number of past bins
%   - n_shuffles: number of surrogates for significance
%
% Outputs:
%   - TE.value: transfer entropy (bits)
%   - TE.p: p-value (shuffle test)
%   - TE.null: null distribution

    if nargin < 4
        n_shuffles = 100;
    end

    assert(length(source) == length(target), 'Inputs must be same length');

    N = length(source);

    % Create history matrices
    X_past = buffer(source(1:end-1), k, k-1, 'nodelay')';
    Y_past = buffer(target(1:end-1), k, k-1, 'nodelay')';
    Y_future = target(k+1:end);

    % Align
    min_len = min([size(X_past,1), size(Y_past,1), length(Y_future)]);

    X_past = X_past(1:min_len, :);
    Y_past = Y_past(1:min_len, :);
    Y_future = Y_future(1:min_len);

    % Build joint states
    joint_past = [Y_past, X_past];

    % Compute TE as conditional mutual information
    TE.value = conditional_mutual_info(Y_future, joint_past, Y_past);

    % Null distribution by shuffling source
    TE.null = zeros(n_shuffles,1);
    for s = 1:n_shuffles
        shuffled = source(randperm(N));
        Xp_s = buffer(shuffled(1:end-1), k, k-1, 'nodelay')';

        % Align lengths
        min_len = min(size(Xp_s,1), size(Y_past,1));
        Xp_s = Xp_s(1:min_len, :);
        Y_past_shuff = Y_past(1:min_len, :);
        Y_future_shuff = Y_future(1:min_len);
        
        % Build surrogate joint past
        jp_s = [Y_past_shuff, Xp_s];
        
        % Compute surrogate TE
        TE.null(s) = conditional_mutual_info(Y_future_shuff, jp_s, Y_past_shuff);
    end

    TE.p = mean(TE.null >= TE.value);  % one-sided p-value
end

function I = conditional_mutual_info(Z, XY, X)
% Compute conditional mutual information I(Z; Y | X)
% Z: future of target (vector)
% XY: [X_past, Y_past]
% X: Y_past only

    z_ids  = vector_to_id(Z);
    xy_ids = matrix_to_id(XY);
    x_ids  = matrix_to_id(X);

    % Probabilities
    p_xyz = joint_hist3(z_ids, xy_ids);
    p_zx  = joint_hist3(z_ids, x_ids);
    p_xy  = accumarray(xy_ids, 1); p_xy = p_xy / sum(p_xy);
    p_x   = accumarray(x_ids, 1);  p_x  = p_x / sum(p_x);

    % Index mapping
    [~,~,b_z]  = unique(z_ids);
    [~,~,b_xy] = unique(xy_ids);
    [~,~,b_x]  = unique(x_ids);

    % Conditional MI
    if any(b_xy > numel(p_xy)) || any(b_x > numel(p_x))
        error('Indexing error: check mapping or histogram computation');
    end
    I = 0;
    for i = 1:length(z_ids)
        p_joint  = p_xyz(b_z(i), b_xy(i));
        p_z_cond = p_zx(b_z(i), b_x(i));
        p_xy_    = p_xy(b_xy(i));
        p_x_     = p_x(b_x(i));

        if p_joint > 0 && p_z_cond > 0 && p_xy_ > 0 && p_x_ > 0
            I = I + p_joint * log2((p_joint * p_x_) / (p_z_cond * p_xy_));
        end
    end
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