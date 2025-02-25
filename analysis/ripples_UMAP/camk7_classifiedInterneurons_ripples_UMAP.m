%% Camk7 sessions (classifying interneurons)
% General configuration
addpath(genpath('C:\Users\pabad\Downloads\umapFileExchange (4.4)\umap'));
cd('C:\Users\pabad\OneDrive - imim.es\Code\structure_index')
SI = py.importlib.import_module('structure_index');
np = py.importlib.import_module('numpy');

n_neighbors = 2;
params = { ...
    'graph_type', 'weighted', ...
    'overlap_method', 'continuity', ...
    'discrete_bin_label', false, ...
    'verbose', true, ...
    'num_shuffles', 10, ...
    'n_bins', 10 ...
};

% Dimensionality reduction parameters
n_neighbors_UMAP = 15;
min_dist = 0.1;
D = 4; % (Intrinsic dimension based on Liset's paper)
metric = 'euclidean';


verbose = true;
sessions = {'Z:\fCamk7\fCamk7_220421_sess17'};
cell_sessions = {'fCamk7_220421_sess17'};

cd('Z:\fCamk7\fCamk7_220421_sess17');
cmap = flip(brewermap([],'RdYlBu'));

% Get interneurons identity

file = dir('*cell_metrics.cellinfo.mat');
load(file.name);

ground_truth = unique(cell_metrics.ground_truth_classification.cell_types);

pv_cells = ismember(cell_metrics.ground_truth_classification.cell_types,'PV+');
sst_cells = ismember(cell_metrics.ground_truth_classification.cell_types,'SST+');
vip_cells = ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+');
id2_cells = ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+');
camk2_cells = ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2');


clear unit

cell_aux = [];
last_cell = 0;
unit = [];

cell_mask = [];

% Compute Ripples UMAP
addpath(genpath('C:\Users\pabad\Downloads\umapFileExchange (4.4)\umap'));

% Dimensionality reduction
n_neighbors = 15;
min_dist = 0.1;
D = 4; % (Intrinsic dimension based on Liset's paper)
metric = 'euclidean';

clear UMAP_ripples
clear spikes_UMAP_v2
clear freq
clear amp


file = dir('*UMAP_ripples_v2.mat');
load(file.name);

UMAP_ripples = ripples_UMAP_v2.maps_filtered;
freq = ripples_UMAP_v2.peakFrequency;
amp = ripples_UMAP_v2.peakAmplitude;

file = dir('*spikes_ripples_UMAP_v2.mat');
load(file.name);

spikes_UMAP_v2 = spikes_ripple_UMAP_v2;

manifold = UMAP_ripples;

try
[reduction,umap,clusterIdentifiers,extras] = run_umap(manifold,...
    'min_dist',min_dist,'n_neighbors',n_neighbors_UMAP,'metric',metric,...
    'n_components',D);

catch
end
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
SI_amp_original= double(amp_SI{1});
SI_amp_shuffled_original = double(amp_SI{4});
SI_amp_original = prctile(SI_amp_shuffled_original,99);
SI_amp_original = prctile(SI_amp_shuffled_original,95);

amp_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(amp), int8(3), params = params );
SI_amp_reduction = double(amp_reduction_SI{1});
SI_amp_reduction_shuffled = double(amp_reduction_SI{4});
SI_amp_reduction = prctile(SI_amp_reduction_shuffled,99);
SI_amp_reduction = prctile(SI_amp_reduction_shuffled,95);


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
SI_freq_original = double(freq_SI{1});
SI_freq_shuffled_original = double(freq_SI{4});
SI_freq_original = prctile(SI_freq_shuffled_original,99);
SI_freq_original = prctile(SI_freq_shuffled_original,95);

freq_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(freq), int8(3), params = params );
SI_freq_reduction = double(freq_reduction_SI{1});
SI_freq_reduction_shuffled = double(freq_reduction_SI{4});
SI_freq_reduction = prctile(SI_freq_reduction_shuffled,99);
SI_freq_reduction = prctile(SI_freq_reduction_shuffled,95);


for ii = 1:length(ground_truth)

    if ~strcmpi(lower(ground_truth{ii}),lower('Noisy_unit'))

        cells = find(ismember(cell_metrics.ground_truth_classification.cell_types,ground_truth{ii}));

        for jj = 1:length(cells)

            y_zscore = zscore(spikes_UMAP_v2(cells(jj),:));

            % Get spikes structure index
            spikes_SI = SI.compute_structure_index(np.array(manifold), np.array(y_zscore), int8(3), params = params );
            SI_spikes_original.(lower(ground_truth{ii}(1:end-1))).values(jj) = double(spikes_SI{1}); 
            SI_spikes_shuffled_original = double(spikes_SI{4});
            SI_spikes_original.(lower(ground_truth{ii}(1:end-1))).prctile99(jj) = prctile(SI_spikes_shuffled_original,99);
            SI_spikes_original.(lower(ground_truth{ii}(1:end-1))).prctile95(jj) = prctile(SI_spikes_shuffled_original,95);
    
            spikes_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(y_zscore), int8(3), params = params );
            SI_spikes_reduction.(lower(ground_truth{ii}(1:end-1))).values(jj) = double(spikes_reduction_SI{1});
            SI_spikes_reduction_shuffled = double(spikes_reduction_SI{4});
            SI_spikes_reduction.(lower(ground_truth{ii}(1:end-1))).prctile99(jj) = prctile(SI_spikes_reduction_shuffled,99);
            SI_spikes_reduction.(lower(ground_truth{ii}(1:end-1))).prctile95(jj) = prctile(SI_spikes_reduction_shuffled,95);

            % disp(SI_spikes_original.(genotype{ii}).(sessions.(genotype{ii}){jj}).value(kk));
            % disp(SI_spikes_reduction.(genotype{ii}){jj}.value(kk));

        end
    end
end

%%

for zz = 1
        cmap = brewermap([],'Spectral');
        colors.camk2 = cmap(5,:);
        colors.camk2_dark =  cmap(5,:)/1.5;
    
        colors.id2 = cmap(100,:);
        colors.id2_dark = cmap(100,:)/1.5;
        
        colors.sst = cmap(200,:);
        colors.sst_dark = cmap(200,:)/1.5;
        
        colors.pv = cmap(240,:);
        colors.pv_dark = cmap(240,:)/2;
        
        cmap = brewermap([],'PiYG');
        colors.vip = cmap(50,:);
        colors.vip_dark = cmap(50,:)/1.5;
        
        colors.vip_ww = cmap(50,:);
        colors.vip_ww_dar = cmap(50,:)/1.5;
    
        cmap = brewermap([],'PRGn');
        colors.vip_nw = cmap(50,:);
        colors.vip_nw_dark = cmap(50,:)/1.5;
        
        cmap = brewermap(100,'RdBu');
        colors.pyr = cmap(25,:);
        colors.pyr_dark = cmap(5,:);
        colors.pyr_light = cmap(40,:);
    
        colors.all5 = [colors.pv; colors.sst; colors.id2; colors.vip; colors.camk2];
        colors.all6 = [colors.pv; colors.sst; colors.id2; colors.vip_ww; colors.vip_nw; colors.camk2];
        
        colors.int = cmap(75,:);
        colors.int_dark = cmap(95,:);
        colors.int_light = cmap(60,:);
end


addpath(genpath('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus'))
figure;

[gs] = groupStats({SI_spikes_original.camk2.values, SI_spikes_original.pv.values, SI_spikes_original.sst.values,...
    SI_spikes_original.vip.values, SI_spikes_original.id2.values},[], ...
    'color',[colors.camk2; colors.pv; colors.sst; colors.vip; colors.id2],'plotData',true,'sigStar',true);
ylim([-0.01 0.15])

[gs] = groupStats({SI_spikes_original.camk2.values(SI_spikes_original.camk2.values > SI_spikes_original.camk2.prctile99), ...
    SI_spikes_original.pv.values(SI_spikes_original.pv.values > SI_spikes_original.pv.prctile99), ...
    SI_spikes_original.sst.values(SI_spikes_original.sst.values > SI_spikes_original.sst.prctile99),...
    SI_spikes_original.vip.values(SI_spikes_original.vip.values > SI_spikes_original.vip.prctile99),...
    SI_spikes_original.id2.values(SI_spikes_original.id2.values > SI_spikes_original.id2.prctile99)},[], ...
    'color',[colors.camk2; colors.pv; colors.sst; colors.vip; colors.id2],'plotData',true,'sigStar',true);
ylim([-0.01 0.15])


[gs] = groupStats({SI_spikes_reduction.camk2.values, SI_spikes_reduction.pv.values, SI_spikes_reduction.sst.values,...
    SI_spikes_reduction.vip.values, SI_spikes_reduction.id2.values},[], ...
    'color',[colors.camk2; colors.pv; colors.sst; colors.vip; colors.id2],'plotData',true,'sigStar',true);
ylim([-0.01 0.15])

cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data\UMAP ripples');
save('camk7_SI_spikes_original_v2.mat','SI_spikes_original');
save('camk7_SI_spikes_reduction_v2.mat','SI_spikes_reduction');

cd('Z:\fCamk7\fCamk7_220421_sess17')
save('camk7_SI_amp_original.mat','SI_amp_original');
save('camk7_SI_amp_reduction.mat','SI_amp_reduction');

save('camk7_SI_freq_original.mat','SI_freq_original');
save('camk7_SI_freq_reduction.mat','SI_freq_reduction');





















