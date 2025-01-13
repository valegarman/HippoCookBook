% General configuration
addpath(genpath('C:\Users\pabad\Downloads\umapFileExchange (4.4)\umap'));
addpath(genpath('C:\Users\pabad\OneDrive - imim.es\Code\HippoCookBook\external_packages\brewerMap_toolbox'));

cd('C:\Users\pabad\OneDrive - imim.es\Code\structure_index')
SI = py.importlib.import_module('structure_index');
np = py.importlib.import_module('numpy');

n_neighbors = 3;
params = { ...
    'graph_type', 'weighted', ...
    'overlap_method', 'continuity', ...
    'discrete_bin_label', false, ...
    'verbose', true, ...
    'num_shuffles', 10, ...
    'n_bins', 10 ...
};

verbose = false;
cmap = flip(brewermap([],'RdYlBu'));

% Dimensionality reduction parameters
n_neighbors_UMAP = 15;
min_dist = 0.1;
D = 4; % (Intrinsic dimension based on Liset's paper)
metric = 'euclidean';

cd('C:\Users\pabad\Desktop')
file = dir('sessions.mat');
load(file.name);

count = 1;
for ii = 1:length(sessions.basepaths)

    cd(sessions.basepaths{ii});

    file = dir('*UMAP_ripples_v2.mat');
    load(file.name);

    UMAP_ripples= ripples_UMAP_v2.maps_filtered;
    freq = ripples_UMAP_v2.peakFrequency;
    amp = ripples_UMAP_v2.peakAmplitude;

    file = dir('*spikes_ripples_UMAP_v2.mat');
    load(file.name);

    spikes_UMAP_v2 = spikes_ripple_UMAP_v2;

    % run UMAP

    manifold = UMAP_ripples;

    try
        [reduction,umap,clusterIdentifiers,extras] = run_umap(manifold,...
        'min_dist',min_dist,'n_neighbors',n_neighbors_UMAP,'metric',metric,...
        'n_components',D);
    catch
        disp('Not possible to run UMAP...');
    end 
    close all;

    % Get the 5th and 95th percentiles of the amplitude data
    vmin = prctile(amp, 5);
    vmax = prctile(amp, 95);

    if verbose
        figure;
        s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
        % Labels and title
        xlabel('UMAP1');
        ylabel('UMAP2');
        title('Amplitude');
        hold on;
        scatter(reduction(:,1), reduction(:,2), 10, amp, 'filled');
        % Set the color scale limits based on the percentiles (similar to vmin and vmax in Python)
        caxis([vmin vmax]);
        colormap(cmap);
        colorbar
    end

    % Get amplitude structure index
    amp_SI = SI.compute_structure_index(np.array(manifold), np.array(amp), int8(3), params = params );
    SI_amp_original.value(ii) = double(amp_SI{1});
    SI_amp_shuffled_original = double(amp_SI{4});
    SI_amp_original.prctile99(ii) = prctile(SI_amp_shuffled_original,99);
    SI_amp_original.prctile95(ii) = prctile(SI_amp_shuffled_original,95);

    amp_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(amp), int8(3), params = params );
    SI_amp_reduction.value(ii) = double(amp_reduction_SI{1});
    SI_amp_reduction_shuffled = double(amp_reduction_SI{4});
    SI_amp_reduction.prctile99(ii) = prctile(SI_amp_reduction_shuffled,99);
    SI_amp_reduction.prctile95(ii) = prctile(SI_amp_reduction_shuffled,95);

    % Get the 5th and 95th percentiles of the frequency data
    vmin = prctile(freq, 5);
    vmax = prctile(freq, 95);

    if verbose
        figure;
        s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
        % Labels and title
        xlabel('UMAP1');
        ylabel('UMAP2');
        title('Frequency');
        hold on;
        scatter(reduction(:,1), reduction(:,2), 10, freq, 'filled');
        % Set the color scale limits based on the percentiles (similar to vmin and vmax in Python)
        caxis([vmin vmax]);
        colormap(cmap);
        colorbar
    end

    % Get frequency structure index
    freq_SI = SI.compute_structure_index(np.array(manifold), np.array(freq), int8(3), params = params );
    SI_freq_original.value(ii) = double(freq_SI{1});
    SI_freq_shuffled_original = double(freq_SI{4});
    SI_freq_original.prctile99(ii) = prctile(SI_freq_shuffled_original,99);
    SI_freq_original.prctile95(ii) = prctile(SI_freq_shuffled_original,95);

    freq_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(freq), int8(3), params = params );
    SI_freq_reduction.value(ii) = double(freq_reduction_SI{1});
    SI_freq_reduction_shuffled = double(freq_reduction_SI{4});
    SI_freq_reduction.prctile99(ii) = prctile(SI_freq_reduction_shuffled,99);
    SI_freq_reduction.prctile95(ii) = prctile(SI_freq_reduction_shuffled,95);

    for jj = 1:size(spikes_UMAP_v2,1)

        y_zscore = zscore(spikes_UMAP_v2(jj,:));
        if verbose
            figure;
            s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
            % Labels and title
            xlabel('UMAP1');
            ylabel('UMAP2');
            title('Spikes');
            hold on;
            scatter(reduction(:,1),reduction(:,2),10,y_zscore,'filled','MarkerFaceAlpha',1);
            caxis([-3 3])
            colormap(cmap);
            colorbar;
        end

        % Get spikes structure index
        spikes_SI = SI.compute_structure_index(np.array(manifold), np.array(y_zscore), int8(3), params = params );
        SI_spikes_original.value(count) = double(spikes_SI{1}); 
        SI_spikes_shuffled_original = double(spikes_SI{4});
        SI_spikes_original.prctile99(count) = prctile(SI_spikes_shuffled_original,99);
        SI_spikes_original.prctile95(count) = prctile(SI_spikes_shuffled_original,95);

        spikes_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(y_zscore), int8(3), params = params );
        SI_spikes_reduction.value(count) = double(spikes_reduction_SI{1});
        SI_spikes_reduction_shuffled = double(spikes_reduction_SI{4});
        SI_spikes_reduction.prctile99(count) = prctile(SI_spikes_reduction_shuffled,99);
        SI_spikes_reduction.prctile95(count) = prctile(SI_spikes_reduction_shuffled,95);

        count = count + 1;
    end
end



cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data');

save('SI_spikes_original.mat','SI_spikes_original');
save('SI_spikes_reduction.mat','SI_spikes_reduction');

save('SI_freq_original.mat','SI_freq_original');
save('SI_freq_reduction.mat','SI_freq_reduction');

save('SI_amp_original.mat','SI_amp_original');
save('SI_amp_reduction.mat','SI_amp_reduction');


