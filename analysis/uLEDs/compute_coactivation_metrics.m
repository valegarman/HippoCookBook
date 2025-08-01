
function coactivation_metrics = compute_coactivation_metrics(uLEDcoactivation, varargin)
% coactivation_metrics = compute_coactivation_metrics(uLEDcoactivation, varargin)
% Compute some metrics from coactivation experiments MV 2025
%
% Parse options
p = inputParser;
addRequired(p,'uLEDcoactivation');
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'save_as','coactivation_metrics',@ischar);
addParameter(p,'include_units', 'all'); % possible values: binary vector, 'pyramidal cells', 'all' (default),
% 'distinct shank', 'distinct_electrode', 'pyramidal cells distinct shank',
% 'pyramidal cells distinct electrode'
addParameter(p,'clipping_threshold',3,@isnumeric);
addParameter(p,'do_plot',true,@islogical);
addParameter(p,'nShuffles',1000,@isnumeric);


parse(p, uLEDcoactivation, varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
force = p.Results.force;
save_as = p.Results.save_as;
include_units = p.Results.include_units;
clipping_threshold = p.Results.clipping_threshold;
do_plot = p.Results.do_plot;
nShuffles = p.Results.nShuffles;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.coactivation_metrics.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Coactivation metrics already computed! Loading file...');
    load(targetFile.name);
    return
end

if ~isfield(uLEDcoactivation, 'clusters') || ~isfield(uLEDcoactivation.clusters, 'artificial_clusters')
    warning('No artificial clusters found; skipping ARI analysis.');
    return
end

if isnumeric(include_units)
    coact = uLEDcoactivation.pairwise.coactZ;
    pre = uLEDcoactivation.pairwise.preZ;
    post = uLEDcoactivation.pairwise.postZ;
    
    coact(~include_units,:) = NaN; coact(:, ~include_units) = NaN;
    pre(~include_units,:) = NaN; pre(:, ~include_units) = NaN;
    post(~include_units,:) = NaN; post(:, ~include_units) = NaN;

elseif ischar(include_units)
    switch lower(include_units)
        case 'all'
            coact = uLEDcoactivation.pairwise.coactZ;
            pre = uLEDcoactivation.pairwise.preZ;
            post = uLEDcoactivation.pairwise.postZ;
        case 'pyramidal cells'
            coact = uLEDcoactivation.pairwise_pyramidalCells.coactZ;
            pre = uLEDcoactivation.pairwise_pyramidalCells.preZ;
            post = uLEDcoactivation.pairwise_pyramidalCells.postZ;
        case 'distinct shank'
            coact = uLEDcoactivation.pairwise_distinct_shank.coactZ;
            pre = uLEDcoactivation.pairwise_distinct_shank.preZ;
            post = uLEDcoactivation.pairwise_distinct_shank.postZ;
        case 'distinct electrode'
            coact = uLEDcoactivation.pairwise_distinct_electrode.coactZ;
            pre = uLEDcoactivation.pairwise_distinct_electrode.preZ;
            post = uLEDcoactivation.pairwise_distinct_electrode.postZ;
        case 'pyramidal cells distinct shank'
            coact = uLEDcoactivation.pairwise_pyramidalCells_distinct_shank.coactZ;
            pre = uLEDcoactivation.pairwise_pyramidalCells_distinct_shank.preZ;
            post = uLEDcoactivation.pairwise_pyramidalCells_distinct_shank.postZ;
        case 'pyramidal cells distinct electrode'
            coact = uLEDcoactivation.pairwise_pyramidalCells_distinct_electrode.coactZ;
            pre = uLEDcoactivation.pairwise_pyramidalCells_distinct_electrode.preZ;
            post = uLEDcoactivation.pairwise_pyramidalCells_distinct_electrode.postZ;
    end
end

%% 1 % CLUSTERING RESPONSES
% 1.1 % Remove bad rows/cols
disp('Performing clustering...');

to_remove = all(isnan(coact), 2);  % 
valid_units = ~to_remove;
Z = coact(valid_units, valid_units);
preZ = pre(valid_units, valid_units);
postZ = post(valid_units, valid_units);

% 1.2 % Clip extreme Z-scores
Z_clipped = max(min(Z, clipping_threshold), -clipping_threshold);
preZ_clipped = max(min(preZ, clipping_threshold), -clipping_threshold);
postZ_clipped = max(min(postZ, clipping_threshold), -clipping_threshold);

% Remove diagonal
Z_clipped(logical(eye(size(Z_clipped)))) = 0;
preZ_clipped(logical(eye(size(preZ_clipped)))) = 0;
postZ_clipped(logical(eye(size(postZ_clipped)))) = 0;

% 1.3 % Signed modularity clustering
% Modularity (Q) measures how well a network is divided into clusters:
% higher Q means stronger within-cluster connectivity vs. random expectation.
% Typical: Q > 0.3 = modular structure; Q > 0.5 = strong community structure.
gamma = 1.1;
method = 'negative_asym';

Q_stim = compare_modularity_to_null(Z_clipped, gamma, nShuffles, method, false);
Ci = Q_stim.clusters_indices;
Q = Q_stim.Q_real;
Q_pre = compare_modularity_to_null(preZ_clipped, gamma, nShuffles, method, false);
Q_post = compare_modularity_to_null(postZ_clipped, gamma, nShuffles, method, false);

Q_pre_by_stim = compare_fixed_modularity_to_null(preZ_clipped, Ci, nShuffles, false);
Q_post_by_stim = compare_fixed_modularity_to_null(postZ_clipped, Ci, nShuffles, false);

fprintf('Modularity (Q) - Pre: %.3f | Stim: %.3f | Post: %.3f\n', Q_pre_by_stim.Q_real, Q, Q_post_by_stim.Q_real);

% 1.4 % Assortability
% Assortativity measures whether nodes preferentially connect to others in the same group.
% Positive = within-group connections; negative = between-group preference; zero = random.
r_pre = compare_assortativity_to_null(preZ_clipped, Ci, nShuffles, false);
r_stim = compare_assortativity_to_null(Z_clipped, Ci, nShuffles, false);
r_post = compare_assortativity_to_null(postZ_clipped, Ci, nShuffles, false);

fprintf('Assortability by cluster - Pre: %.3f | Stim: %.3f | Post: %.3f\n', r_pre.r_real, r_stim.r_real, r_post.r_real);

% save values
coactivation_metrics.community_louvain.Q_pre = Q_pre;
coactivation_metrics.community_louvain.Q_stim = Q_stim;
coactivation_metrics.community_louvain.Q_post = Q_post;

coactivation_metrics.community_louvain.Q_pre_by_stim = Q_pre_by_stim;
coactivation_metrics.community_louvain.Q_post_by_stim = Q_post_by_stim;

coactivation_metrics.community_louvain.cluster_stim = Ci;
coactivation_metrics.community_louvain.cluster_pre = Q_pre.clusters_indices;
coactivation_metrics.community_louvain.cluster_post = Q_post.clusters_indices;

coactivation_metrics.assortability.r_pre = r_pre;
coactivation_metrics.assortability.r_stim = r_stim;
coactivation_metrics.assortability.r_post = r_post;

if do_plot
    % stim clusters
    caxis_lim = 2;
    figure
    subplot(3,3,5)
    [~, perm] = sort(Ci);
    sorted_matrix = Z_clipped(perm, perm);
    sorted_Ci = Ci(perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    title({['Stim'], ['Q = ' num2str(Q,2) ', r = ' num2str(r_stim.r_real,2)]});
    hold off;
    axis square
    
    subplot(3,3,4)
    sorted_matrix = preZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Pre'], ['Q = ' num2str(Q_pre_by_stim.Q_real,2) ', r = ' num2str(r_pre.r_real,2)]});
    ylabel('Clustering by Stim', 'FontWeight','bold');
    axis square
    
    subplot(3,3,6)
    sorted_matrix = postZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Post'], ['Q = ' num2str(Q_post_by_stim.Q_real,2) ', r = ' num2str(r_post.r_real,2)]});
    axis square
    
    % pre clusters
    cluster_pre = Q_pre.clusters_indices;
    caxis_lim = 2;
    subplot(3,3,2)
    [~, perm] = sort(cluster_pre);
    sorted_matrix = Z_clipped(perm, perm);
    sorted_Ci = cluster_pre(perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Stim']});
    axis square
    
    subplot(3,3,1)
    sorted_matrix = preZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Pre'], ['Q = ' num2str(Q_pre.Q_real,2)]});
    ylabel('Clustering by Pre', 'FontWeight','bold');
    axis square
    
    subplot(3,3,3)
    sorted_matrix = postZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Post']});
    axis square

    % post clusters
    cluster_post = Q_post.clusters_indices;
    caxis_lim = 2;
    subplot(3,3,8)
    [~, perm] = sort(cluster_post);
    sorted_matrix = Z_clipped(perm, perm);
    sorted_Ci = cluster_post(perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Stim']});
    xlabel('# Neurons');
    axis square
    
    subplot(3,3,7)
    sorted_matrix = preZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Pre']});
    ylabel('Clustering by Post', 'FontWeight','bold');
    axis square
    
    subplot(3,3,9)
    sorted_matrix = postZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Post'], ['Q = ' num2str(Q_post.Q_real,2)]});
    axis square

    mkdir('SummaryFigures');
    exportgraphics(gcf,['SummaryFigures\' save_as '_louvain_clusters.png']);
end


%% 2 % COMPARING PARTITIONS BY RAND INDEX
% Adjusted Rand Index (ARI) quantifies similarity between two clusterings.
% ARI = 1: perfect match; ARI = 0: random agreement; ARI < 0: worse than chance.
% Used here to compare artificial labels with Louvain clustering.
disp('Comparing particions...');

artificial_clusters = uLEDcoactivation.clusters.artificial_clusters;
validIdx = ~isnan(artificial_clusters) & artificial_clusters ~= 0;
true_labels  = artificial_clusters(validIdx);
pre_labels   = Q_pre.clusters_indices(validIdx);
stim_labels  = Q_stim.clusters_indices(validIdx);
post_labels  = Q_post.clusters_indices(validIdx);

ari_pre  = adjusted_rand_index(true_labels, pre_labels);
ari_stim = adjusted_rand_index(true_labels, stim_labels);
ari_post = adjusted_rand_index(true_labels, post_labels);
ari_pre_post = adjusted_rand_index(pre_labels, post_labels);
ari_pre_stim = adjusted_rand_index(pre_labels, stim_labels);
ari_post_stim = adjusted_rand_index(post_labels, stim_labels);

fprintf('ARI - Pre: %.3f | Stim: %.3f | Post: %.3f\n', ari_pre, ari_stim, ari_post);
fprintf('ARI - Pre-Post: %.3f | Pre-Stim: %.3f | Post-Stim: %.3f\n', ari_pre_post, ari_pre_stim, ari_post_stim);

nShuffles2 = 10000;
ari_shuff_stim = zeros(nShuffles2,1);
ari_shuff_pre = zeros(nShuffles2,1);
ari_shuff_post = zeros(nShuffles2,1);
ari_shuff_pre_post = zeros(nShuffles2,1);
ari_shuff_pre_stim = zeros(nShuffles2,1);
ari_shuff_post_stim = zeros(nShuffles2,1);

for i = 1:nShuffles2
    shuffled = true_labels(randperm(length(true_labels)));
    ari_shuff_stim(i) = adjusted_rand_index(shuffled, stim_labels);
    ari_shuff_pre(i) = adjusted_rand_index(shuffled, pre_labels);
    ari_shuff_post(i) = adjusted_rand_index(shuffled, post_labels);
    %
    ari_shuff_pre_post(i) = adjusted_rand_index(shuffled, post_labels);
    ari_shuff_pre_stim(i) = adjusted_rand_index(shuffled, pre_labels);
    ari_shuff_post_stim(i) = adjusted_rand_index(shuffled, post_labels);
end

real_ari = [ari_pre, ari_stim, ari_post, ari_pre_post, ari_pre_stim, ari_post_stim];
shuffled_ari = {ari_shuff_pre, ari_shuff_stim, ari_shuff_post, ari_shuff_pre_post, ari_shuff_pre_stim, ari_shuff_post_stim};
labels = {'Pre', 'Stim', 'Post', 'Pre-Post', 'Pre-Stim', 'Stim-Post'};

% Zscore and p-value
z_scores = zeros(1,3);
p_vals = zeros(1,3);

for i = 1:6
    mu = mean(shuffled_ari{i});
    sigma = std(shuffled_ari{i});
    z_scores(i) = (real_ari(i) - mu) / sigma;
    p_vals(i) = mean(shuffled_ari{i} >= real_ari(i));  % unilateral (efecto positivo)
end

% save values
coactivation_metrics.artificial_clusters = artificial_clusters';
coactivation_metrics.adjusted_rand_index.ari_pre = ari_pre;
coactivation_metrics.adjusted_rand_index.ari_stim = ari_stim;
coactivation_metrics.adjusted_rand_index.ari_post = ari_post;
coactivation_metrics.adjusted_rand_index.ari_pre_post = ari_pre_post;
coactivation_metrics.adjusted_rand_index.ari_pre_stim = ari_pre_stim;
coactivation_metrics.adjusted_rand_index.ari_stim_post = ari_post_stim;

coactivation_metrics.adjusted_rand_index.ari_pre_pval = p_vals(1);
coactivation_metrics.adjusted_rand_index.ari_stim_pval = p_vals(2);
coactivation_metrics.adjusted_rand_index.ari_post_pval = p_vals(3);
coactivation_metrics.adjusted_rand_index.ari_pre_post_pval = p_vals(4);
coactivation_metrics.adjusted_rand_index.ari_pre_stim_pval = p_vals(5);
coactivation_metrics.adjusted_rand_index.ari_stim_post_pval = p_vals(6);

coactivation_metrics.adjusted_rand_index.ari_pre_Z = z_scores(1);
coactivation_metrics.adjusted_rand_index.ari_stim_Z = z_scores(2);
coactivation_metrics.adjusted_rand_index.ari_post_Z = z_scores(3);
coactivation_metrics.adjusted_rand_index.ari_pre_post_Z = z_scores(4);
coactivation_metrics.adjusted_rand_index.ari_pre_stim_Z = z_scores(5);
coactivation_metrics.adjusted_rand_index.ari_stim_post_Z = z_scores(6);

if do_plot

    figure('Position', [300 300 800 600]);
    all_vals = cat(1, shuffled_ari{:}, real_ari(:));
    xlims = [floor(min(all_vals)*100)/100, ceil(max(all_vals)*100)/100];
    for i = 1:3
        subplot(3,3,i)
        histogram(shuffled_ari{i}, 40, 'FaceColor', [0.8 0.8 0.8], ...
                  'EdgeColor', 'none', 'Normalization', 'count');
        hold on
        xline(real_ari(i), 'r-', 'LineWidth', 2);
        xlim(xlims)
        title([labels{i}, '  —  Z = ', num2str(z_scores(i), '%.2f'), ...
               ',  p = ', num2str(p_vals(i), '%.3f')])
        if i == 1
            ylabel('Count')
        end
    end
    
    subplot(3,3,4:6)
    bar(z_scores, 'FaceColor', [0.5 0.7 0.9])
    set(gca, 'XTickLabel', labels, 'FontSize', 10)
    ylabel('Z-score (vs. shuffle)')
    ylim([0 max(z_scores)*1.2])
    title('Z-scores of Observed ARI')
    for i = 1:3
        text(i, z_scores(i)+0.05, ['p = ', num2str(p_vals(i), '%.3f')], ...
            'HorizontalAlignment', 'center', 'FontSize', 10)
    end
    sgtitle('Observed ARI vs. Shuffling');

    caxis_lim = 2;
    subplot(3,3,8)
    [~, perm] = sort(artificial_clusters);
    sorted_matrix = Z_clipped(perm, perm);
    sorted_Ci = artificial_clusters(perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Stim']});
    xlabel('# Neurons');
    axis square

    subplot(3,3,7)
    sorted_matrix = preZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Pre']});
    ylabel('Clustering by uLEDs groups', 'FontWeight','bold');
    axis square
    
    subplot(3,3,9)
    sorted_matrix = postZ_clipped(perm, perm);
    imagesc(sorted_matrix); colormap('parula'); colorbar;
    caxis([-caxis_lim caxis_lim]);
    hold on;
    b = find(diff(sorted_Ci) ~= 0) + 0.5;
    for x = b; xline(x,'k','LineWidth',1.5); yline(x,'k','LineWidth',1.5); end
    hold off;
    title({['Post']});
    axis square

    exportgraphics(gcf,['SummaryFigures\' save_as '_adjusted_rand_index.png']);

    % Sankey diagrams
    output_name = ['sankey_comparison'];

    % Pre vs Artificial
    T1 = table(Q_pre.clusters_indices(:), artificial_clusters(:), 'VariableNames', {'Ci1', 'Ci2'});
    writetable(T1, 'sankey_pre.csv');
    
    % Stim vs Artificial
    T2 = table(Q_stim.clusters_indices(:), artificial_clusters(:), 'VariableNames', {'Ci1', 'Ci2'});
    writetable(T2, 'sankey_stim.csv');
    
    % Post vs Artificial
    T3 = table(Q_pre.clusters_indices(:), artificial_clusters(:), 'VariableNames', {'Ci1', 'Ci2'});
    writetable(T3, 'sankey_post.csv');
    
    script_path = which('sankey_triple_plot.py');
    python_path = char(pyenv().Executable);
    command = ['"', python_path, '" "', script_path, '" ', output_name];
    system(command);

    img_path = fullfile('SummaryFigures', [output_name, '.png']);
    img = imread(img_path);
    fig = figure('Name', ['Sankey: ', output_name]);
    imshow(img);

    % Sankey multi step, pre, stim, post
    output_name = 'sankey_pre_stim_post';

    T = table(Q_pre.clusters_indices(:), Q_stim.clusters_indices(:), Q_post.clusters_indices(:), ...
              'VariableNames', {'Pre', 'Stim', 'Post'});
    writetable(T, 'multi_stage_clusters.csv');
    
    script_path = which('sankey_multistage_plot.py');
    python_path = char(pyenv().Executable);
    command = ['"', python_path, '" "', script_path, '" ', output_name];
    system(command);
    
    img_path = fullfile('SummaryFigures', [output_name, '.png']);
    img = imread(img_path);
    fig = figure('Name', ['Sankey: ', output_name]);
    imshow(img);

    % Sankey multi step, pre, artificial, post
    output_name = 'sankey_artificial_stim_post';

    T = table(Q_pre.clusters_indices(:), artificial_clusters(:), Q_post.clusters_indices(:), ...
              'VariableNames', {'Pre', 'Artificial', 'Post'});
    writetable(T, 'multi_stage_clusters.csv');
    
    script_path = which('sankey_multistage_plot.py');
    python_path = char(pyenv().Executable);
    command = ['"', python_path, '" "', script_path, '" ', output_name];
    system(command);
    
    img_path = fullfile('SummaryFigures', [output_name, '.png']);
    img = imread(img_path);
    fig = figure('Name', ['Sankey: ', output_name]);
    imshow(img);
end

%% 3. ESTIMATING NEURON CONTRIBUTIONS
disp('Estimating neuron contributions to particions...');

% 3.1 % Leave-One-Out ARI
C_all = [Q_pre.clusters_indices(:), Q_stim.clusters_indices(:), Q_post.clusters_indices(:)];
C_ref = artificial_clusters(:);

loo_ari = compute_loo_ari_all(C_all, C_ref);

% 3.2 % Pairwise consistency by cell pair and familiy
pairwise_consistency = compute_pairwise_consistency_all(C_all, C_ref);

[cell_types, cell_classification_stats, cell_metrics, cell_subtypes] = cellTypeClassifier('modelType','hippocampus5');

% compute change of consistency by cell family
A_pre  = pairwise_consistency.full_matrices{1};  % Pre
A_post = pairwise_consistency.full_matrices{3};  % Post
delta_A = A_post - A_pre;  % -1: lost, +1: won, 0: no change
n = size(delta_A,1);

[unique_types, ~, type_idx] = unique(cell_types);  % type_idx = [1, 2, 3, ...]
n_types = numel(unique_types);

delta_sum = zeros(n_types,1);
n_pairs   = zeros(n_types,1);

for i = 1:n
    type_i = type_idx(i);  % cell class of neuron i
    for j = 1:n
        if i == j, continue; end
        delta = delta_A(i,j);  % -1, 0, o 1
        delta_sum(type_i) = delta_sum(type_i) + delta;
        n_pairs(type_i) = n_pairs(type_i) + 1;
    end
end

delta_mean = delta_sum ./ n_pairs;  % consistency

delta_shuff = zeros(n_types, nShuffles);

for s = 1:nShuffles
    shuffled_idx = type_idx(randperm(length(type_idx)));  % shuffling

    delta_sum = zeros(n_types,1);
    n_pairs = zeros(n_types,1);

    for i = 1:n
        type_i = shuffled_idx(i);
        for j = 1:n
            if i == j, continue; end
            delta = delta_A(i,j);
            delta_sum(type_i) = delta_sum(type_i) + delta;
            n_pairs(type_i) = n_pairs(type_i) + 1;
        end
    end

    delta_shuff(:,s) = delta_sum ./ n_pairs;
end

p_vals = zeros(n_types,1);
for t = 1:n_types
    real_val = delta_mean(t);
    null_dist = delta_shuff(t,:);
    p_vals(t) = mean(abs(null_dist) >= abs(real_val));  % dos colas
    z_score_val(t) = (real_val - mean(null_dist)) / std(null_dist);
end

pairwise_consistency.consistency_post_pre_between_families.z_score_over_null = z_score_val;
pairwise_consistency.consistency_post_pre_between_families.p_vals = p_vals;
pairwise_consistency.consistency_post_pre_between_families.post_pre_normalized = delta_mean;
pairwise_consistency.consistency_post_pre_between_families.post_pre = delta_sum;
pairwise_consistency.consistency_post_pre_between_families.n_pairs = n_pairs;

% save values
coactivation_metrics.adjusted_rand_index.leave_one_out_test = loo_ari;
coactivation_metrics.pairwise_consistency = pairwise_consistency;

if do_plot
    % Unique classes (excluding 'Noisy_unit')
    cat_types = categorical(cell_types);
    [~, ~, numeric_labels] = unique(cat_types);
    type_names = categories(cat_types);
    colors = getColors(type_names);
    
    % importance
    importance = loo_ari.importance_zscore;
    figure
    subplot(3,3,1)
    gs = groupStats(importance(:,2), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Stim','FontWeight','normal');
    ylabel('Leave-one-out importance (SD)');
    
    subplot(3,3,2)
    gs = groupStats(importance(:,3), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post','FontWeight','normal');
    
    subplot(3,3,3)
    gs = groupStats(importance(:,3) - importance(:,1), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post-pre','FontWeight','normal');
    
    % pairwise consistency
    consistency = pairwise_consistency.zscore;
    subplot(3,3,4)
    gs = groupStats(consistency(:,2), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Stim','FontWeight','normal');
    ylabel('Pairwise consistency (SD)');
    
    subplot(3,3,5)
    gs = groupStats(consistency(:,3), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post','FontWeight','normal');
    
    subplot(3,3,6)
    gs = groupStats(consistency(:,3) - consistency(:,1), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post-pre','FontWeight','normal');
    
    % within vs between cluster consistency
    within_between = pairwise_consistency.within_vs_between_zscore;
    subplot(3,3,7)
    gs = groupStats(within_between(:,2), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Stim','FontWeight','normal');
    ylabel('Within vs between cluster consistency (SD)');
    
    subplot(3,3,8)
    gs = groupStats(within_between(:,3), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post','FontWeight','normal');
    
    subplot(3,3,9)
    gs = groupStats(within_between(:,3) - within_between(:,1), numeric_labels, 'plotType', 'roundPlot', 'color', colors, 'plotData', true, 'inAxis', true);
    set(gca, 'XTick', 1:length(type_names), 'XTickLabel', type_names, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title('Post-pre','FontWeight','normal');

    % subplot(3,3,7)
    % b = bar(z_score_val);
    % b.FaceColor = 'flat';
    % b.CData = getColors(unique_types);
    % set(gca, 'XTick', 1:n_types, 'XTickLabel', unique_types, ...
    %          'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none')
    % ylabel('Δ Pairwise consistency against null (post - pre)')
    % hold on
    % sig_idx = find(p_vals < 0.05);
    % for i = sig_idx'
    %     text(i, delta_mean(i) + 0.01, '*', 'FontSize', 16, ...
    %          'HorizontalAlignment', 'center', 'Color', 'k')
    % end

    exportgraphics(gcf,['SummaryFigures\' save_as '_clustering_neuron_contribution.png']);
end

%% 4. CLUSTER LEVEL METRICS
disp('Computing cluster-wise metrics...');

% 4.1 By stim clusters
[cluster_interactions_stim] = computeClusterInteractions(Z_clipped, Q_stim.clusters_indices);
[cluster_interactions_pre_by_stim] = computeClusterInteractions(preZ_clipped, Q_stim.clusters_indices);
[cluster_interactions_post_by_stim] = computeClusterInteractions(postZ_clipped, Q_stim.clusters_indices);

% Signed LOO impact: negative values suggest the neuron functionally promoted excitation or coordination across clusters.
% Positive values suggest the neuron suppressed or inhibited interactions between clusters.
[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(Z_clipped, Q_stim.clusters_indices, nShuffles, false);
loo_interactions_stim = computeLOOimpact_allMetrics_withNull(Z_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(preZ_clipped, Q_stim.clusters_indices, nShuffles, false);
loo_interactions_pre_by_stim = computeLOOimpact_allMetrics_withNull(preZ_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(postZ_clipped, Q_stim.clusters_indices, nShuffles, false);
loo_interactions_post_by_stim = computeLOOimpact_allMetrics_withNull(postZ_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

loo_interactions_post_pre_by_stim = subtractLOOstructures(loo_interactions_post_by_stim, loo_interactions_pre_by_stim);

% save values
coactivation_metrics.cluster_metrics_by_stim.cluster_interactions_stim = cluster_interactions_stim;
coactivation_metrics.cluster_metrics_by_stim.cluster_interactions_pre_by_stim = cluster_interactions_pre_by_stim;
coactivation_metrics.cluster_metrics_by_stim.cluster_interactions_post_by_stim = cluster_interactions_post_by_stim;

coactivation_metrics.cluster_metrics_by_stim.loo_interactions_stim = loo_interactions_stim;
coactivation_metrics.cluster_metrics_by_stim.loo_interactions_pre_by_stim = loo_interactions_pre_by_stim;
coactivation_metrics.cluster_metrics_by_stim.loo_interactions_post_by_stim = loo_interactions_post_by_stim;
coactivation_metrics.cluster_metrics_by_stim.loo_interactions_post_pre_by_stim = loo_interactions_post_pre_by_stim;

% 4.2 By artificial clusters
unique_vals = unique(artificial_clusters(~isnan(artificial_clusters)));  % si hay NaN los ignora
Ci_mapped = zeros(size(artificial_clusters));
for k = 1:length(unique_vals)
    Ci_mapped(artificial_clusters == unique_vals(k)) = k;
end
artificial_clusters2 = Ci_mapped;
[cluster_interactions_stim_by_artificial] = computeClusterInteractions(Z_clipped, artificial_clusters2);
[cluster_interactions_pre_by_artificial] = computeClusterInteractions(preZ_clipped, artificial_clusters2);
[cluster_interactions_post_by_artificial] = computeClusterInteractions(postZ_clipped, artificial_clusters2);

% Signed LOO impact: negative values suggest the neuron functionally promoted excitation or coordination across clusters.
% Positive values suggest the neuron suppressed or inhibited interactions between clusters.
[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(Z_clipped, artificial_clusters2, nShuffles, false);
loo_interactions_stim_by_artificial = computeLOOimpact_allMetrics_withNull(Z_clipped, artificial_clusters2, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(preZ_clipped,artificial_clusters2, nShuffles, false);
loo_interactions_pre_by_artificial = computeLOOimpact_allMetrics_withNull(preZ_clipped, artificial_clusters2, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(postZ_clipped, artificial_clusters2, nShuffles, false);
loo_interactions_post_by_artificial = computeLOOimpact_allMetrics_withNull(postZ_clipped, artificial_clusters2, Z_nulls, Ci_nulls);

loo_interactions_post_pre_by_artificial = subtractLOOstructures(loo_interactions_post_by_stim, loo_interactions_pre_by_stim);

% save values
coactivation_metrics.cluster_metrics_by_artificial.cluster_interactions_stim_by_artificial = cluster_interactions_stim_by_artificial;
coactivation_metrics.cluster_metrics_by_artificial.cluster_interactions_pre_by_artificial = cluster_interactions_pre_by_artificial;
coactivation_metrics.cluster_metrics_by_artificial.cluster_interactions_post_by_artificial = cluster_interactions_post_by_artificial;

coactivation_metrics.cluster_metrics_by_artificial.loo_interactions_stim_by_artificial = loo_interactions_stim_by_artificial;
coactivation_metrics.cluster_metrics_by_artificial.loo_interactions_pre_by_artificial = loo_interactions_pre_by_artificial;
coactivation_metrics.cluster_metrics_by_artificial.loo_interactions_post_by_artificial = loo_interactions_post_by_artificial;
coactivation_metrics.cluster_metrics_by_artificial.loo_interactions_post_pre_by_artificial = loo_interactions_post_pre_by_artificial;

if do_plot
    % plots by stim
    figure
    subplot(2,3,1)
    gs = groupStats({cluster_interactions_pre_by_stim.diagValues(:), cluster_interactions_stim.diagValues(:),...
        cluster_interactions_post_by_stim.diagValues(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Mean Interaction by stim(SD)');
    title('Intra-cluster','FontWeight','normal');
    
    subplot(2,3,2)
    gs = groupStats({cluster_interactions_pre_by_stim.neighborPairs(:), cluster_interactions_stim.neighborPairs(:),...
        cluster_interactions_post_by_stim.neighborPairs(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    title('Neigh-cluster','FontWeight','normal');
    
    subplot(2,3,3)
    gs = groupStats({cluster_interactions_pre_by_stim.offDiag(:), cluster_interactions_stim.offDiag(:),...
        cluster_interactions_post_by_stim.offDiag(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    title('Inter-cluster','FontWeight','normal');


    % plots by artificial clusters
    subplot(2,3,4)
    gs = groupStats({cluster_interactions_pre_by_artificial.diagValues(:), cluster_interactions_stim_by_artificial.diagValues(:),...
        cluster_interactions_post_by_artificial.diagValues(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Mean Interaction by artificial (SD)');
    title('Intra-cluster','FontWeight','normal');
    
    subplot(2,3,5)
    gs = groupStats({cluster_interactions_pre_by_artificial.neighborPairs(:), cluster_interactions_stim_by_artificial.neighborPairs(:),...
        cluster_interactions_post_by_artificial.neighborPairs(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    title('Neigh-cluster','FontWeight','normal');
    
    subplot(2,3,6)
    gs = groupStats({cluster_interactions_pre_by_artificial.offDiag(:), cluster_interactions_stim_by_artificial.offDiag(:),...
        cluster_interactions_post_by_artificial.offDiag(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    title('Inter-cluster','FontWeight','normal');

    exportgraphics(gcf,['SummaryFigures\' save_as '_cluster_level_metrics.png']);
end

%% 5. NEURON-WISE CLUSTER METRICS
disp('Computing neuron-wise metrics...');

% intra: Median interaction of each neuron with other neurons from its own cluster
%        → High values indicate strong coupling within its group (functional coherence).
%
% inter: Median interaction with neurons from all other clusters (excluding own)
%        → High values suggest cross-cluster communication or influence.
%
% neigh: Median interaction with neurons from neighboring clusters (e.g., cluster ±1)
%        → Useful if clusters have a meaningful ordering (e.g., cortical layers).
%
% participation_total: Normalized entropy of the absolute interaction vector across all clusters
%        → Measures how broadly a neuron interacts across the network (regardless of sign).
%        → High = uniform influence across clusters; low = focused on a few.
%
% participation_excitation: Normalized entropy of the positive interaction vector
%        → Measures how broadly a neuron distributes its excitatory influence.
%        → High = excites many clusters; low = excites mainly one.
%
% participation_inhibition: Difference between total and excitation participation
%        → Captures the diversity of inhibitory (negative) interactions.
%        → High = broad inhibitory influence; low = focused inhibition or none.
% 4.1 By stim clusters
[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(preZ_clipped, Q_stim.clusters_indices, nShuffles, false);
cluster_cell_metrics_pre_by_stim = computeNeuronClusterMetricsWithNull(preZ_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(Z_clipped, Q_stim.clusters_indices, nShuffles, false);
cluster_cell_metrics_stim = computeNeuronClusterMetricsWithNull(Z_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(postZ_clipped, Q_stim.clusters_indices, nShuffles, false);
cluster_cell_metrics_post_by_stim = computeNeuronClusterMetricsWithNull(postZ_clipped, Q_stim.clusters_indices, Z_nulls, Ci_nulls);

cluster_cell_metrics_post_pre_by_stim = subtractMetricsStructures(cluster_cell_metrics_post_by_stim, cluster_cell_metrics_pre_by_stim);

% save values
coactivation_metrics.neuron_metrics_by_stim.cluster_cell_metrics_pre_by_stim = cluster_cell_metrics_pre_by_stim;
coactivation_metrics.neuron_metrics_by_stim.cluster_cell_metrics_stim = cluster_cell_metrics_stim;
coactivation_metrics.neuron_metrics_by_stim.cluster_cell_metrics_post_by_stim = cluster_cell_metrics_post_by_stim;
coactivation_metrics.neuron_metrics_by_stim.cluster_cell_metrics_post_pre_by_stim = cluster_cell_metrics_post_pre_by_stim;

% 4.1 By artificial clusters
[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(preZ_clipped, artificial_clusters2', nShuffles, false);
cluster_cell_metrics_pre_by_artificial = computeNeuronClusterMetricsWithNull(preZ_clipped, artificial_clusters2', Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(Z_clipped, artificial_clusters2', nShuffles, false);
cluster_cell_metrics_stim_by_artificial = computeNeuronClusterMetricsWithNull(Z_clipped, artificial_clusters2', Z_nulls, Ci_nulls);

[Z_nulls, Ci_nulls] = generateClusterPermutationNulls(postZ_clipped, artificial_clusters2', nShuffles, false);
cluster_cell_metrics_post_by_artificial = computeNeuronClusterMetricsWithNull(postZ_clipped, artificial_clusters2', Z_nulls, Ci_nulls);

cluster_cell_metrics_post_pre_by_artificial = subtractMetricsStructures(cluster_cell_metrics_post_by_artificial, cluster_cell_metrics_pre_by_artificial);

% save values
coactivation_metrics.neuron_metrics_by_artificial.cluster_cell_metrics_pre_by_artificial = cluster_cell_metrics_pre_by_artificial;
coactivation_metrics.neuron_metrics_by_artificial.cluster_cell_metrics_stim_by_artificial = cluster_cell_metrics_stim_by_artificial;
coactivation_metrics.neuron_metrics_by_artificial.cluster_cell_metrics_post_by_artificial = cluster_cell_metrics_post_by_artificial;
coactivation_metrics.neuron_metrics_by_artificial.cluster_cell_metrics_post_pre_by_artificial = cluster_cell_metrics_post_pre_by_artificial;

% 4.2 By stim cluster, separating possitive and negative interactions
Z_pos = max(Z_clipped, 0);
Z_neg = min(Z_clipped, 0);

if do_plot
    % plots by stim
    figure
    subplot(2,3,1)
    gs = groupStats({cluster_cell_metrics_pre_by_stim.z_score.intra(:), cluster_cell_metrics_stim.z_score.intra(:),...
        cluster_cell_metrics_post_by_stim.z_score.intra(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Median interaction by stim (SD)');
    title('Intra-cluster','FontWeight','normal');

    subplot(2,3,2)
    gs = groupStats({cluster_cell_metrics_pre_by_stim.z_score.neigh(:), cluster_cell_metrics_stim.z_score.neigh(:),...
        cluster_cell_metrics_post_by_stim.z_score.neigh(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Median interaction by stim (SD)');
    title('Neigh-cluster','FontWeight','normal');

    subplot(2,3,3)
    gs = groupStats({cluster_cell_metrics_post_pre_by_stim.z_score.intra(:),...
        cluster_cell_metrics_post_pre_by_stim.z_score.neigh(:), cluster_cell_metrics_post_pre_by_stim.z_score.participation_total(:), ...
        cluster_cell_metrics_post_pre_by_stim.z_score.participation_excitation(:), cluster_cell_metrics_post_pre_by_stim.z_score.participation_inhibition(:)}, [], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', false);
    ylabel('\DeltaMedian interaction by stim (SD)');
    ylim([-5 5]);
    set(gca, 'XTick', [1:5], 'XTickLabel', {'Intra', 'Neigh', 'Participation', 'Part Exc', 'Part Inh'});

    subplot(2,3,4)
    gs = groupStats({cluster_cell_metrics_pre_by_artificial.z_score.intra(:), cluster_cell_metrics_stim_by_artificial.z_score.intra(:),...
        cluster_cell_metrics_post_by_artificial.z_score.intra(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Median interaction by stim (SD)');
    title('Intra-cluster','FontWeight','normal');

    subplot(2,3,5)
    gs = groupStats({cluster_cell_metrics_pre_by_artificial.z_score.neigh(:), cluster_cell_metrics_stim_by_artificial.z_score.neigh(:),...
        cluster_cell_metrics_post_by_artificial.z_score.neigh(:)},[], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', true);
    ylabel('Median interaction by stim (SD)');
    title('Neigh-cluster','FontWeight','normal');

    subplot(2,3,6)
    gs = groupStats({cluster_cell_metrics_post_pre_by_artificial.z_score.intra(:),...
        cluster_cell_metrics_post_pre_by_artificial.z_score.neigh(:), cluster_cell_metrics_post_pre_by_artificial.z_score.participation_total(:), ...
        cluster_cell_metrics_post_pre_by_artificial.z_score.participation_excitation(:), cluster_cell_metrics_post_pre_by_artificial.z_score.participation_inhibition(:)}, [], 'plotType', 'roundPlot', 'inAxis', true, 'repeatedMeasures', true, 'plotData', true, 'plotConnectors', false);
    ylabel('\DeltaMedian interaction by stim (SD)');
    ylim([-5 5]);
    set(gca, 'XTick', [1:5], 'XTickLabel', {'Intra', 'Neigh', 'Participation', 'Part Exc', 'Part Inh'});

end

% save results...
if saveMat
    save(fullfile(basepath, [save_as '.coactivation_metrics.cellinfo.mat']), 'coactivation_metrics');
end


cd(prevPath);
end

function result = compute_loo_ari_all(C_all, C_ref)
% Compute LOO-ARI importance per neuron for multiple conditions.
%
% INPUTS:
%   - C_all: N x K matrix of clustering labels (e.g., pre, stim, post)
%   - C_ref: N x 1 vector of reference labels (e.g., artificial clustering)
%
% OUTPUT:
%   - result: struct with fields:
%       * importance: N x K raw LOO-ARI values
%       * importance_relative: N x K (relative to full ARI)
%       * importance_zscore: N x K z-scored values (per condition)
%       * ari_full: 1 x K vector of ARI values using all neurons
%
% Requires: rand_index_adjusted(C1, C2)

nNeurons = size(C_all, 1);
nConds = size(C_all, 2);

importance = zeros(nNeurons, nConds);
importance_relative = zeros(nNeurons, nConds);
importance_zscore = zeros(nNeurons, nConds);
ari_full = zeros(1, nConds);

for k = 1:nConds
    Ci1 = C_all(:,k);
    Ci2 = C_ref(:);
    
    ari_full(k) = adjusted_rand_index(Ci1, Ci2);
    
    for i = 1:nNeurons
        keep = true(nNeurons, 1);
        keep(i) = false;
        ari_reduced = adjusted_rand_index(Ci1(keep), Ci2(keep));
        importance(i,k) = ari_full(k) - ari_reduced;
    end
    
    % Relative importance (avoiding division by 0)
    if ari_full(k) ~= 0
        importance_relative(:,k) = importance(:,k) / ari_full(k);
    else
        importance_relative(:,k) = zeros(nNeurons,1);
    end
    
    % Z-score normalization
    mu = mean(importance(:,k));
    sigma = std(importance(:,k));
    if sigma ~= 0
        importance_zscore(:,k) = (importance(:,k) - mu) / sigma;
    else
        importance_zscore(:,k) = zeros(nNeurons,1);
    end
end

% Output struct
result.importance = importance;
result.importance_relative = importance_relative;
result.importance_zscore = importance_zscore;
result.ari_full = ari_full;

end

function pairwise_consistency = compute_pairwise_consistency_all(C_all, C_ref, nShuffles)
% Compute pairwise consistency and structure metrics vs. a reference clustering
%
% INPUTS:
%   - C_all: N x K matrix of clustering labels (e.g., Pre, Stim, Post)
%   - C_ref: N x 1 reference clustering
%   - nShuffles: optional (default: 1000) number of shuffles for null distributions
%
% OUTPUTS: struct with fields:
%   - consistency: N x K normalized pairwise consistency [0–1]
%   - raw:         N x K raw agreement counts
%   - zscore:      N x K z-scores vs shuffled null
%   - pval:        N x K empirical p-values
%   - full_matrices: 1 x K cell of N x N agreement matrices
%   - within_cluster, between_cluster: N x K mean consistency with same/different clusters
%   - within_vs_between_diff: N x K difference (within - between)
%   - within_vs_between_zscore, within_vs_between_pval: N x K stats from shuffled null

if nargin < 3
    nShuffles = 1000;
end

nNeurons = size(C_all, 1);
nConds = size(C_all, 2);

% Outputs
consistency_raw     = zeros(nNeurons, nConds);
consistency_norm    = zeros(nNeurons, nConds);
consistency_zscore  = zeros(nNeurons, nConds);
p_vals              = zeros(nNeurons, nConds);
agreement_matrices  = cell(1, nConds);

within_cluster      = zeros(nNeurons, nConds);
between_cluster     = zeros(nNeurons, nConds);
within_vs_between_diff = zeros(nNeurons, nConds);
within_vs_between_zscore = zeros(nNeurons, nConds);
within_vs_between_pval   = zeros(nNeurons, nConds);

shuffle_vals = zeros(nNeurons, nConds, nShuffles);
shuffle_diff = zeros(nNeurons, nConds, nShuffles);

for k = 1:nConds
    Ci1 = C_all(:,k);
    Ci2 = C_ref(:); % fixed reference

    A = zeros(nNeurons);
    for i = 1:nNeurons
        for j = i+1:nNeurons
            same1 = Ci1(i) == Ci1(j);
            same2 = Ci2(i) == Ci2(j);
            agree = double(same1 == same2);
            A(i,j) = agree;
            A(j,i) = agree;
        end
    end
    agreement_matrices{k} = A;

    % Real values
    raw = sum(A, 2);
    norm = raw / (nNeurons - 1);

    consistency_raw(:,k)  = raw;
    consistency_norm(:,k) = norm;

    for i = 1:nNeurons
        same_mask = (Ci2 == Ci2(i));
        same_mask(i) = false;
        diff_mask = ~same_mask;

        within_cluster(i,k)  = mean(A(i,same_mask));
        between_cluster(i,k) = mean(A(i,diff_mask));
        within_vs_between_diff(i,k) = within_cluster(i,k) - between_cluster(i,k);
    end

    % Shuffling
    for s = 1:nShuffles
        Ci2_shuf = Ci2(randperm(nNeurons));
        A_shuf = zeros(nNeurons);

        for i = 1:nNeurons
            for j = i+1:nNeurons
                same1 = Ci1(i) == Ci1(j);
                same2 = Ci2_shuf(i) == Ci2_shuf(j);
                agree = double(same1 == same2);
                A_shuf(i,j) = agree;
                A_shuf(j,i) = agree;
            end
        end

        shuffle_vals(:,k,s) = sum(A_shuf, 2) / (nNeurons - 1);

        for i = 1:nNeurons
            same_mask = (Ci2_shuf == Ci2_shuf(i));
            same_mask(i) = false;
            diff_mask = ~same_mask;

            w = mean(A_shuf(i,same_mask));
            b = mean(A_shuf(i,diff_mask));
            shuffle_diff(i,k,s) = w - b;
        end
    end

    % Z-scores and p-values
    null_mean = mean(shuffle_vals(:,k,:), 3);
    null_std  = std(shuffle_vals(:,k,:), 0, 3);
    consistency_zscore(:,k) = (consistency_norm(:,k) - null_mean) ./ null_std;

    for i = 1:nNeurons
        p_vals(i,k) = (1 + sum(shuffle_vals(i,k,:) >= consistency_norm(i,k))) / (1 + nShuffles);

        % Within-between shuffle stats
        real_diff = within_vs_between_diff(i,k);
        null_diff = squeeze(shuffle_diff(i,k,:));
        mu_null = mean(null_diff);
        sigma_null = std(null_diff);
        within_vs_between_zscore(i,k) = (real_diff - mu_null) / sigma_null;
        within_vs_between_pval(i,k) = (1 + sum(null_diff >= real_diff)) / (1 + nShuffles);
    end
end

% Output struct
pairwise_consistency.consistency    = consistency_norm;
pairwise_consistency.raw           = consistency_raw;
pairwise_consistency.zscore        = consistency_zscore;
pairwise_consistency.pval          = p_vals;
pairwise_consistency.full_matrices = agreement_matrices;

pairwise_consistency.within_cluster  = within_cluster;
pairwise_consistency.between_cluster = between_cluster;
pairwise_consistency.within_vs_between_diff    = within_vs_between_diff;
pairwise_consistency.within_vs_between_zscore  = within_vs_between_zscore;
pairwise_consistency.within_vs_between_pval    = within_vs_between_pval;

end

function result = compare_modularity_to_null(Z, gamma, nShuffles, method, doPlot)
% Compare modularity Q to a null distribution.
%
% Inputs:
%   - Z         : NxN matrix (e.g., Z_clipped), can be signed
%   - gamma     : resolution parameter for Louvain (e.g., 1.1)
%   - nShuffles : number of null samples (e.g., 1000)
%   - method    : modularity method ('negative_asym', etc.)
%   - doPlot    : (optional) if true, plot null dist vs Q
%
% Output:
%   - result: struct with fields:
%       .Q_real       : modularity of original Z
%       .Q_null       : vector of shuffled Qs
%       .z_score      : z-score of Q_real vs null
%       .p_value      : p-value (one-tailed)
%       .Q_mean_null  : mean of null Q
%       .Q_std_null   : std of null Q

if nargin < 5
    doPlot = false;
end

% Compute real modularity
rng(1);  % reproducible
[ci, Q_real] = community_louvain(Z, gamma, [], method);

% Preallocate null distribution
Q_null = nan(nShuffles, 1);
N = size(Z,1);

for s = 1:nShuffles
    perm = randperm(N);
    Z_shuffled = Z(perm, perm);
    [~, Q_null(s)] = community_louvain(Z_shuffled, gamma, [], method);
end

% Stats
Q_mean = mean(Q_null);
Q_std = std(Q_null);
z = (Q_real - Q_mean) / Q_std;
p = mean(Q_null >= Q_real);  % one-tailed

% Store results
result.Q_real = Q_real;
result.Q_null = Q_null;
result.z_score = z;
result.p_value = p;
result.Q_mean_null = Q_mean;
result.Q_std_null = Q_std;
result.clusters_indices = ci;


% Optional plot
if doPlot
    histogram(Q_null, 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
    hold on
    y = ylim;
    plot([Q_real Q_real], y, 'r', 'LineWidth', 2);
    hold off
    xlabel('Modularity Q (null)');
    ylabel('Frequency');
    title(sprintf('Q_{real} = %.3f | z = %.2f | p = %.4f', Q_real, z, p));
end
end

function result = compare_fixed_modularity_to_null(Z, Ci, nShuffles, doPlot)
% Compare fixed modularity (Q = modularity_signed(Z, Ci)) to a null.
%
% Inputs:
%   - Z         : NxN connectivity matrix (e.g. preZ_clipped)
%   - Ci        : 1xN or Nx1 vector of community assignments (from stim)
%   - nShuffles : number of null samples (e.g., 1000)
%   - doPlot    : (optional) true/false for plotting null dist
%
% Output:
%   - result: struct with fields:
%       .Q_real       : modularity_signed(Z, Ci)
%       .Q_null       : Q from shuffled matrices (Ci fixed)
%       .z_score      : z-score of Q_real vs null
%       .p_value      : p-value (one-tailed)
%       .Q_mean_null  : mean of null Q
%       .Q_std_null   : std of null Q

if nargin < 4
    doPlot = false;
end

% Ensure column vector
Ci = Ci(:);
N = size(Z,1);

% Real Q
Q_real = modularity_signed(Z, Ci);

% Null distribution
Q_null = nan(nShuffles, 1);
for s = 1:nShuffles
    perm = randperm(N);
    Z_shuffled = Z(perm, perm);  % preserve symmetry structure
    Q_null(s) = modularity_signed(Z_shuffled, Ci);
end

% Stats
Q_mean = mean(Q_null);
Q_std = std(Q_null);
z = (Q_real - Q_mean) / Q_std;
p = mean(Q_null >= Q_real);  % one-tailed test

% Output
result.Q_real = Q_real;
result.Q_null = Q_null;
result.z_score = z;
result.p_value = p;
result.Q_mean_null = Q_mean;
result.Q_std_null = Q_std;

% Plot
if doPlot
    histogram(Q_null, 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
    hold on
    y = ylim;
    plot([Q_real Q_real], y, 'r', 'LineWidth', 2);
    hold off
    xlabel('Modularity Q (null)');
    ylabel('Frequency');
    title(sprintf('Q_{real} = %.3f | z = %.2f | p = %.4f', Q_real, z, p));
end
end

function r = assortativity_by_cluster(W, Ci)
    n = length(Ci);
    same = Ci(:) == Ci(:)';
    upper = triu(true(n), 1);  % upper triangle to avoid double-counting
    weights = W(upper);
    same_clust = same(upper);
    
    r = corr(double(same_clust), weights, 'Type', 'Spearman');  % or Pearson
end

function result = compare_assortativity_to_null(W, Ci, nShuffles, doPlot)
% Compare assortativity_by_cluster(W, Ci) to a null model.
% Shuffles W, keeps Ci fixed.
%
% Inputs:
%   - W         : NxN symmetric matrix of weights (can be signed)
%   - Ci        : community assignments (Nx1 or 1xN)
%   - nShuffles : number of null samples (e.g. 1000)
%   - doPlot    : (optional) if true, plots null distribution vs real
%
% Output:
%   - result: struct with fields:
%       .r_real      : observed assortativity
%       .r_null      : vector of null assortativities
%       .z_score     : z-score of r_real vs null
%       .p_value     : one-tailed p-value (real > null)
%       .r_mean_null : mean of null
%       .r_std_null  : std of null

if nargin < 4
    doPlot = false;
end

% Compute real assortativity
r_real = assortativity_by_cluster(W, Ci);

% Null distribution: shuffle W
N = size(W,1);
r_null = nan(nShuffles,1);

for s = 1:nShuffles
    % Shuffle upper triangle
    upper = triu(true(N), 1);
    w_vals = W(upper);
    w_shuffled = w_vals(randperm(length(w_vals)));
    
    % Rebuild symmetric matrix
    W_shuff = zeros(N);
    W_shuff(upper) = w_shuffled;
    W_shuff = W_shuff + W_shuff';  % symmetric

    % Compute assortativity
    r_null(s) = assortativity_by_cluster(W_shuff, Ci);
end

% Statistics
r_mean = mean(r_null);
r_std = std(r_null);
z = (r_real - r_mean) / r_std;
p = mean(r_null >= r_real);  % one-tailed

% Output
result.r_real = r_real;
result.r_null = r_null;
result.z_score = z;
result.p_value = p;
result.r_mean_null = r_mean;
result.r_std_null = r_std;

% Optional plot
if doPlot
    histogram(r_null, 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
    hold on
    y = ylim;
    plot([r_real r_real], y, 'r', 'LineWidth', 2);
    hold off
    xlabel('Assortativity (null)');
    ylabel('Frequency');
    title(sprintf('r_{real} = %.3f | z = %.2f | p = %.4f', r_real, z, p));
end
end

function [cluster_interactions] = computeClusterInteractions(Z, Ci)
    nClusters = max(Ci);
    interactionMatrix = zeros(nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            cells_i = find(Ci == i);
            cells_j = find(Ci == j);
            submatrix = Z(cells_i, cells_j);
            if i == j
                mask = ~eye(length(cells_i));
                submatrix = submatrix(mask);
            end
            interactionMatrix(i,j) = median(submatrix(:), 'omitnan');
        end
    end

    % Diagonal (intra-cluster)
    diagValues = diag(interactionMatrix);

    % Off-diagonal (inter-cluster)
    temp = interactionMatrix';
    temp(1:nClusters+1:end) = NaN;
    offDiag = temp(~isnan(temp));

    % Between neigh clusters 
    neighborPairs = [];
    for c = 1:nClusters-1
        neighborPairs = [neighborPairs, interactionMatrix(c,c+1), interactionMatrix(c+1,c)];
    end
    neighborPairs = neighborPairs';

    cluster_interactions.interactionMatrix = interactionMatrix;
    cluster_interactions.diagValues = diagValues;
    cluster_interactions.offDiag = offDiag;
    cluster_interactions.neighborPairs = neighborPairs;
end

function loo_results = computeLOOimpact_allMetrics(Z, Ci)
% Computes LOO impact of each neuron on all cluster interaction metrics:
% diagValues, offDiag, neighborPairs
%
% INPUT:
%   - Z: NxN interaction matrix
%   - Ci: Nx1 cluster assignment vector
%
% OUTPUT:
%   - loo_results: struct with fields:
%       .signed_diag (Nx1): signed change in diag values
%       .signed_offDiag (Nx1): signed change in offDiag
%       .signed_neighbor (Nx1): signed change in neighborPairs
%       .loo_diag (Nx1): L2 norm of diag change (abs score)
%       .loo_offDiag (Nx1): L2 norm of offDiag change
%       .loo_neighbor (Nx1): L2 norm of neighborPairs change
%       .valid (Nx1): true if neuron removal was valid

    nNeurons = length(Ci);
    loo_results.signed_diag = nan(nNeurons, 1);
    loo_results.signed_offDiag = nan(nNeurons, 1);
    loo_results.signed_neighbor = nan(nNeurons, 1);

    loo_results.loo_diag = nan(nNeurons, 1);
    loo_results.loo_offDiag = nan(nNeurons, 1);
    loo_results.loo_neighbor = nan(nNeurons, 1);

    loo_results.valid = false(nNeurons, 1);

    % Matriz de referencia con todas las neuronas
    ref = computeClusterInteractions(Z, Ci);
    ref_diag = ref.diagValues;
    ref_offDiag = ref.offDiag;
    ref_neigh = ref.neighborPairs;

    for n = 1:nNeurons
        keep = true(nNeurons, 1);
        keep(n) = false;

        Z_wo = Z(keep, keep);
        Ci_wo = Ci(keep);

        try
            res = computeClusterInteractions(Z_wo, Ci_wo);

            % Validación: debe tener misma dimensión que referencia
            if length(res.diagValues) ~= length(ref_diag)
                continue
            end

            % ------------------
            % DIAGONAL
            delta_diag = res.diagValues - ref_diag;
            loo_results.signed_diag(n) = mean(delta_diag, 'omitnan');
            loo_results.loo_diag(n) = norm(delta_diag, 2);

            % OFF-DIAGONAL
            delta_off = res.offDiag - ref_offDiag;
            loo_results.signed_offDiag(n) = mean(delta_off, 'omitnan');
            loo_results.loo_offDiag(n) = norm(delta_off, 2);

            % NEIGHBORS
            delta_neigh = res.neighborPairs - ref_neigh;
            loo_results.signed_neighbor(n) = mean(delta_neigh, 'omitnan');
            loo_results.loo_neighbor(n) = norm(delta_neigh, 2);

            loo_results.valid(n) = true;

        catch
            continue  % se queda como NaN si falla
        end
    end
end

function loo_results = computeLOOimpact_allMetrics_withNull(Z, Ci, Z_nulls, Ci_nulls)
% Computes LOO impact + null distribution comparison for each neuron and metric
% Includes both absolute (loo) and signed metrics
%
% Inputs:
%   - Z: NxN real interaction matrix
%   - Ci: Nx1 real cluster labels
%   - Z_nulls: cell array of null Z matrices (1 x nShuffles)
%   - Ci_nulls: cell array of null cluster labels (1 x nShuffles)
%
% Output:
%   - loo_results: struct with fields:
%       .signed_*, .loo_*, .z_loo_*, .p_loo_*, .z_signed_*, .p_signed_* (Nx1)
%       .valid_* (Nx1)

    nNeurons = length(Ci);
    nShuffles = length(Z_nulls);
    metrics = {'diag', 'offDiag', 'neighbor'};

    % Real LOO computation
    loo_obs = computeLOOimpact_allMetrics(Z, Ci);

    % Loop over metrics
    for m = 1:length(metrics)
        metric = metrics{m};

        % Observed values
        signed_obs = loo_obs.(['signed_' metric]);
        loo_obs_metric = loo_obs.(['loo_' metric]);
        valid = loo_obs.valid;

        loo_results.(['signed_' metric]) = signed_obs;
        loo_results.(['loo_' metric]) = loo_obs_metric;
        loo_results.(['valid_' metric]) = valid;

        % Allocate null distributions
        null_signed = nan(nNeurons, nShuffles);
        null_loo = nan(nNeurons, nShuffles);

        for s = 1:nShuffles
            Zs = Z_nulls{s};
            Cis = Ci_nulls{s};
            loo_null = computeLOOimpact_allMetrics(Zs, Cis);
            if isfield(loo_null, ['loo_' metric])
                null_signed(:,s) = loo_null.(['signed_' metric]);
                null_loo(:,s) = loo_null.(['loo_' metric]);
            end
        end

        % Z-score and p-value for absolute LOO
        mu_loo = mean(null_loo, 2, 'omitnan');
        std_loo = std(null_loo, 0, 2, 'omitnan');
        loo_results.(['z_loo_' metric]) = (loo_obs_metric - mu_loo) ./ std_loo;

        count_loo = sum(null_loo >= loo_obs_metric, 2, 'omitnan');
        loo_results.(['p_loo_' metric]) = (count_loo + 1) ./ (nShuffles + 1);

        % Z-score and p-value for signed impact
        mu_signed = mean(null_signed, 2, 'omitnan');
        std_signed = std(null_signed, 0, 2, 'omitnan');
        loo_results.(['z_signed_' metric]) = (signed_obs - mu_signed) ./ std_signed;

        count_signed = sum(null_signed >= signed_obs, 2, 'omitnan');
        loo_results.(['p_signed_' metric]) = (count_signed + 1) ./ (nShuffles + 1);
    end
end

function plotLOO_lollipop(loo_score, cell_type)
    % Lollipop plot: each neuron has a vertical line ending in a colored circle (by cell type)
    % Inputs:
    %   - loo_score: Nx1 vector (e.g., LOO impact)
    %   - cell_type: Nx1 categorical, string, or numeric

    loo_score = loo_score(:); % column
    N = length(loo_score);

    % Colors
    [G, ~, group_idx] = unique(cell_type);
    cmap = getColors(G);

    % Make plot
    hold on
    for i = 1:N
        y = loo_score(i);
        if isnan(y), continue; end
        c = cmap(group_idx(i), :);
        line([i i], [0 y], 'Color', c, 'LineWidth', 1.2);
        plot(i, y, 'o', 'MarkerSize', 6, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end

    xlim([0 N+1]);
    ylabel('LOO Impact');
    xlabel('Neuron index');
    title('Neuron-wise LOO impact (colored by cell type)');

    % Legend
    legend_entries = arrayfun(@(g) sprintf('%s', string(G(g))), 1:length(G), 'UniformOutput', false);
    dummy = arrayfun(@(k) plot(nan, nan, 'o', 'Color', cmap(k,:), 'MarkerFaceColor', cmap(k,:)), 1:length(G));
    legend(dummy, legend_entries, 'Location', 'northeastoutside');
end

function delta_struct = subtractLOOstructures(post_struct, pre_struct)
% Subtracts all matching fields from two LOO result structures
% delta_struct = post_struct - pre_struct (field by field)

    delta_struct = struct();

    fields = fieldnames(post_struct);

    for i = 1:length(fields)
        f = fields{i};
        if isnumeric(post_struct.(f)) && isnumeric(pre_struct.(f)) ...
                && isequal(size(post_struct.(f)), size(pre_struct.(f)))
            delta_struct.(f) = post_struct.(f) - pre_struct.(f);
        else
            delta_struct.(f) = [];  % leave empty if not compatible
        end
    end
end

function metrics = computeNeuronClusterMetrics(Z, Ci)
% Computes neuron-wise cluster interaction metrics from a signed interaction matrix Z and cluster labels Ci.
% Includes intra/inter/neighbor interaction medians and entropy-based participation indices.
% Internally computes positive and negative matrices (Z_pos, Z_neg) for excitation/inhibition separation.

nCells = length(Ci);
nClusters = max(Ci);

% Prepare positive (excitation-like) and negative (inhibition-like) interaction matrices
Z_pos = max(Z, 0);       % positive values only
Z_neg = abs(min(Z, 0));  % absolute of negative values

% Initialize metrics
intra_mean = NaN(nCells, 1);
inter_mean = NaN(nCells, 1);
neighbor_mean = NaN(nCells, 1);
p_matrix_pos = zeros(nCells, nClusters);
p_matrix_neg = zeros(nCells, nClusters);

for i = 1:nCells
    my_cluster = Ci(i);
    same_cluster = find(Ci == my_cluster & (1:nCells)' ~= i);
    other_cluster = find(Ci ~= my_cluster);

    % Neighboring cluster indices
    neighbors = [];
    if my_cluster > 1
        neighbors = [neighbors, my_cluster - 1];
    end
    if my_cluster < nClusters
        neighbors = [neighbors, my_cluster + 1];
    end
    idx_neighbors = find(ismember(Ci, neighbors));

    % Metric calculations
    if ~isempty(same_cluster)
        intra_mean(i) = median(Z(i, same_cluster), 'omitnan');
    end
    if ~isempty(other_cluster)
        inter_mean(i) = median(Z(i, other_cluster), 'omitnan');
    end
    if ~isempty(idx_neighbors)
        neighbor_mean(i) = median(Z(i, idx_neighbors), 'omitnan');
    end

    for c = 1:nClusters
        idx_c = find(Ci == c);
        p_matrix_pos(i, c) = mean(Z_pos(i, idx_c), 'omitnan');
        p_matrix_neg(i, c) = mean(Z_neg(i, idx_c), 'omitnan');
    end
end

% Compute entropy-based participation indices

% Excitation
row_sums_pos = sum(p_matrix_pos, 2);
p_norm_pos = p_matrix_pos ./ (row_sums_pos + eps);
entropy_pos = -sum(p_norm_pos .* log(p_norm_pos + eps), 2);
participation_excitation = entropy_pos / log(nClusters);

% Inhibition
row_sums_neg = sum(p_matrix_neg, 2);
p_norm_neg = p_matrix_neg ./ (row_sums_neg + eps);
entropy_neg = -sum(p_norm_neg .* log(p_norm_neg + eps), 2);
participation_inhibition = entropy_neg / log(nClusters);

% Total
p_matrix_total = p_matrix_pos + p_matrix_neg;
row_sums_total = sum(p_matrix_total, 2);
p_norm_total = p_matrix_total ./ (row_sums_total + eps);
entropy_total = -sum(p_norm_total .* log(p_norm_total + eps), 2);
participation_total = entropy_total / log(nClusters);

% Output
metrics.intra = intra_mean;
metrics.inter = inter_mean;
metrics.neigh = neighbor_mean;
metrics.participation_excitation = participation_excitation;
metrics.participation_inhibition = participation_inhibition;
metrics.participation_total = participation_total;
metrics.Z_pos = Z_pos;
metrics.Z_neg = Z_neg;
end

function metrics = computeNeuronClusterMetricsWithNull(Z, Ci, Z_nulls, Ci_nulls)
% Computes neuron-wise cluster metrics + null distribution z-scores and p-values.
% This version internally supports metrics derived from both excitation and inhibition.
%
% Inputs:
%   - Z: NxN interaction matrix (with signs)
%   - Ci: Nx1 cluster labels
%   - Z_nulls: {1 x S} shuffled Z matrices
%   - Ci_nulls: {1 x S} shuffled cluster labels
%
% Output:
%   - metrics.real.(field)
%   - metrics.null.(field)
%   - metrics.z_score.(field)
%   - metrics.p_value.(field)

    % Compute observed metrics for full, positive, and negative interactions
    metrics_obs_full = computeNeuronClusterMetrics(Z, Ci);
    metrics_obs_pos = computeNeuronClusterMetrics(max(Z, 0), Ci);
    metrics_obs_neg = computeNeuronClusterMetrics(abs(min(Z, 0)), Ci);

    % Define metric fields (exclude raw matrices)
    metric_fields = setdiff(fieldnames(metrics_obs_full), {'Z_pos', 'Z_neg'});

    types = {'', '_pos', '_neg'};
    obs_structs = {metrics_obs_full, metrics_obs_pos, metrics_obs_neg};

    % Store real values for each type
    for t = 1:numel(types)
        suffix = types{t};
        obs = obs_structs{t};
        for f = 1:numel(metric_fields)
            field = metric_fields{f};
            metrics.real.([field suffix]) = obs.(field);
        end
    end

    nNeurons = length(Ci);
    nShuffles = length(Z_nulls);

    % Null distributions, z-scores, p-values for each metric and type
    for f = 1:numel(metric_fields)
        field = metric_fields{f};

        % Initialize
        null_values = nan(nNeurons, nShuffles);
        null_values_pos = nan(nNeurons, nShuffles);
        null_values_neg = nan(nNeurons, nShuffles);

        for s = 1:nShuffles
            Zs = Z_nulls{s};
            Cis = Ci_nulls{s};
            m_full = computeNeuronClusterMetrics(Zs, Cis);
            m_pos = computeNeuronClusterMetrics(max(Zs, 0), Cis);
            m_neg = computeNeuronClusterMetrics(abs(min(Zs, 0)), Cis);

            null_values(:, s)     = m_full.(field);
            null_values_pos(:, s) = m_pos.(field);
            null_values_neg(:, s) = m_neg.(field);
        end

        % Store nulls
        metrics.null.(field)     = null_values;
        metrics.null.([field '_pos']) = null_values_pos;
        metrics.null.([field '_neg']) = null_values_neg;

        % Compute z-score and p-value
        for suffix = { '', '_pos', '_neg' }
            sfx = suffix{1};
            nulls = metrics.null.([field sfx]);
            real_vals = metrics.real.([field sfx]);

            mu = mean(nulls, 2, 'omitnan');
            sigma = std(nulls, 0, 2, 'omitnan');

            metrics.z_score.([field sfx]) = (real_vals - mu) ./ (sigma + eps);
            metrics.p_value.([field sfx]) = (sum(nulls >= real_vals, 2) + 1) ./ (nShuffles + 1);
        end
    end
end

function [p_total, p_exc, p_inh, P_abs, P_pos] = computeParticipationIndices(Z, Ci, nClusters)
% Returns participation indices and raw interaction matrices per neuron

    nCells = length(Ci);
    P = zeros(nCells, nClusters);  % neuron × cluster

    for c = 1:nClusters
        idx = find(Ci == c);
        P(:, c) = mean(Z(:, idx), 2);
    end

    % Matrix versions
    P_abs = abs(P);
    P_pos = max(P, 0);

    % Entropies
    p_total = entropyIndex(P_abs);
    p_exc   = entropyIndex(P_pos);
    p_inh   = p_total - p_exc;
    p_inh(p_inh == 0) = NaN;
end

function H = entropyIndex(M)
% Computes row-wise normalized Shannon entropy
    row_sums = sum(M, 2);
    P = M ./ (row_sums + eps);
    H = -sum(P .* log(P + eps), 2) / log(size(M, 2));
end

function [Z_nulls, Ci_nulls] = generateClusterPermutationNulls(Z, Ci, nShuffles, permute_Z)
% Generate shuffled null distributions of Ci (and optionally Z)
%
% Inputs:
%   - Z: NxN interaction matrix
%   - Ci: Nx1 cluster labels
%   - nShuffles: number of permutations
%   - permute_Z: logical, whether to permute Z rows/cols (default: false)
%
% Outputs:
%   - Z_nulls: {1 x nShuffles} cell array of NxN matrices
%   - Ci_nulls: {1 x nShuffles} cell array of Nx1 vectors

    if nargin < 4
        permute_Z = false;
    end

    N = length(Ci);
    Z_nulls = cell(1, nShuffles);
    Ci_nulls = cell(1, nShuffles);

    for s = 1:nShuffles
        perm = randperm(N);
        Ci_nulls{s} = Ci(perm);
        if permute_Z
            Z_nulls{s} = Z(perm, perm);
        else
            Z_nulls{s} = Z;
        end
    end
end

function delta = subtractMetricsStructures(metrics_post, metrics_pre)
% Subtracts all matching subfields between two full metric structures (e.g., 'real', 'z_score', etc.)
% Returns a structure delta with same hierarchy: delta.real.intra, etc.
% Skips fields that are not column vectors (e.g. Z_pos, null, etc.)

    top_fields = intersect(fieldnames(metrics_post), fieldnames(metrics_pre));
    delta = struct();

    for i = 1:numel(top_fields)
        tf = top_fields{i};
        if isstruct(metrics_post.(tf)) && isstruct(metrics_pre.(tf))
            subfields = intersect(fieldnames(metrics_post.(tf)), fieldnames(metrics_pre.(tf)));
            for j = 1:numel(subfields)
                sf = subfields{j};
                a = metrics_post.(tf).(sf);
                b = metrics_pre.(tf).(sf);

                % Only subtract if both are vectors of same size
                if isvector(a) && isvector(b) && isequal(size(a), size(b))
                    delta.(tf).(sf) = a - b;
                end
            end
        end
    end
end