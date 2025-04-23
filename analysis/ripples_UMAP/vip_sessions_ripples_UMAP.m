% cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data')
% file = dir('chapter_0_data19-Sep-2024.mat');
% load(file.name);

saveFig = false;
%% VIP SESSIONS
cmap = flip(brewermap([],'RdYlBu'));
sessions = {'Z:\fCr1\fCr1_220301_sess3','Z:\fCr1\fCr1_220304_sess6','Z:\fCr1\fCr1_220322_sess17'};
% sessions = {'Z:\fCr1\fCr1_220304_sess6'};

clear unit
% sst_cell = find(taggedCells.hippo.sst(ismember(projectResults.session,lower({'fSst1_190604_sess34','fSst1_190605_sess35','fSst1_190606_sess36'}))));

% cell_sessions = {'fCr1_220304_sess6'};
cell_sessions = {'fCr1_220301_sess3','fCr1_220304_sess6','fCr1_220322_sess17'};

cell_aux = [];
last_cell = 0;
unit = [];

for ii = 1:length(cell_sessions)
    cell_aux = find(taggedCells.hippo.vip(ismember(projectResults.session,lower(cell_sessions{ii}))));
    cell_aux = cell_aux + last_cell;
    unit = [unit cell_aux];
    last_cell = last_cell + projectSessionResults.numcells(find(ismember(lower(projectSessionResults.sessionName),lower(cell_sessions{ii}))));
end

cell_mask = [];
for ii = 1:length(sessions)
    sess = strsplit(sessions{ii},'\');
    sess = sess{3};
    % sst_cell{ii} = find(taggedCells.hippo.sst(ismember(projectResults.session,lower(sst_sess))));

    cell_mask = [cell_mask ones(1,length(find(taggedCells.hippo.vip(ismember(projectResults.session,lower(sess))))))*ii];
end

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

mask = [];
for ii = 1:length(sessions)

    cd(sessions{ii})
    
    file = dir('*UMAP_ripples_v2.mat');
    load(file.name);

    UMAP_ripples{ii} = ripples_UMAP_v2.maps_filtered;
    freq{ii} = ripples_UMAP_v2.peakFrequency;
    amp{ii} = ripples_UMAP_v2.peakAmplitude;
    
    file = dir('*spikes_ripples_UMAP_v2.mat');
    load(file.name);

    spikes_UMAP_v2{ii} = spikes_ripple_UMAP_v2;

    % Mask for ripples
   
    mask = [mask ones(1,length(freq{ii}))*ii];

end


manifold = cell2mat(UMAP_ripples');

try
[reduction,umap,clusterIdentifiers,extras] = run_umap(manifold,...
    'min_dist',min_dist,'n_neighbors',n_neighbors,'metric',metric,...
    'n_components',D);
catch
end
close all;

%% AMPLITUDE PROJECTED
% Assuming 'ripple_data.amp' is the amplitude data used for coloring
amp = cell2mat(amp');

% Get the 5th and 95th percentiles of the amplitude data
vmin = prctile(amp, 5);
vmax = prctile(amp, 95);

% Create the scatter plot with the colors based on the amplitude data
% Project frequency and amplitude
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

% Set the colormap to 'coolwarm' equivalent in MATLAB. MATLAB doesn't have a 'coolwarm' colormap,
% but you can use 'cool' or 'parula', or download a custom 'coolwarm' colormap.
colormap(cmap);  % 'cool' is similar to 'coolwarm' in matplotlib
colorbar

if saveFig
    cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP');
    saveas(s,'vip_amp_ripples_UMAP.fig');
end

%% FREQUENCY PROJECTED

% Assuming 'ripple_data.amp' is the amplitude data used for coloring
freq = cell2mat(freq');

% Get the 5th and 95th percentiles of the amplitude data
vmin = prctile(freq, 5);
vmax = prctile(freq, 95);

% Create the scatter plot with the colors based on the amplitude data
% Project frequency and amplitude
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

% Set the colormap to 'coolwarm' equivalent in MATLAB. MATLAB doesn't have a 'coolwarm' colormap,
% but you can use 'cool' or 'parula', or download a custom 'coolwarm' colormap.
colormap(cmap);  % 'cool' is similar to 'coolwarm' in matplotlib
colorbar

if saveFig
    cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP');
    saveas(s,'vip_freq_ripples_UMAP.fig');
end

%% SPIKES
clear y
for ii = 1:length(spikes_UMAP_v2)
    y{ii} = zscore(spikes_UMAP_v2{ii},[],2);
end

a = max(cell2mat(cellfun(@length,y,'UniformOutput',false)));

for ii = 1:length(y)
    if size(y{ii},2) ~= a
        y{ii}(:,end:a) = NaN;
    end
end

y_out = cell2mat(y');
y_zscore = y_out;
color_indices = ceil(rescale(y_zscore, 1, 100));

vmin = prctile(y_zscore(unit,:), 5,'all');
vmax = prctile(y_zscore(unit,:), 95,'all');

cmin = -max(max(y_zscore(unit,:)));
cmax = max(max(y_zscore(unit,:)));

cmap = flip(brewermap([],'RdYlBu'));

% Let's plot separated
for ii = 1:length(unit)

    figure;
    s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;        
    scatter(reduction(mask==cell_mask(ii),1),reduction(mask==cell_mask(ii),2),10,y_zscore(unit(ii),1:1:length(find(mask==cell_mask(ii)))),'filled','MarkerFaceAlpha',.8);
    colormap(cmap);
    % caxis([0 5])
    % caxis([min(min(y_zscore(sst_cell,:))) max(max(y_zscore(sst_cell,:)))])
    % caxis([vmin vmax])
    % caxis([-1 1])
    caxis([-.5 .5])
    colorbar;
end

%% Compute the mean for each session when more than one neuron
a = max(cell2mat(cellfun(@length,spikes_UMAP_v2,'UniformOutput',false)));

for ii = 1:length(spikes_UMAP_v2)
    if size(spikes_UMAP_v2{ii},2) ~= a
        spikes_UMAP_v2{ii}(:,end:a) = NaN;
    end
end

y = cell2mat(spikes_UMAP_v2');
y_mean_zscored = nan(length(unique(cell_mask)),size(y,2));
y_mean = nan(length(unique(cell_mask)),size(y,2));

for ii = 1:length(unique(cell_mask))
    if length(unit(cell_mask == ii)) == 1
        y_mean(ii,:) = y(unit(cell_mask == ii),:);
    else
        y_mean(ii,:) = mean(y(unit(cell_mask == ii),:));
    end

    y_mean_zscored(ii,1:length(find(~isnan(y_mean(ii,:))))) = zscore(y_mean(ii,~isnan(y_mean(ii,:))));
    if any(isnan(y_mean(ii,:)))
        y_mean_zscored(ii,length(find(~isnan(y_mean(ii,:)))):end) = nan;
    end
    % y_mean_zscored(length(find(~isnan(y_mean(ii,:)))):end) = nan;
end


vmin = prctile(y_mean_zscored, 5,'all');
vmax = prctile(y_mean_zscored, 95,'all');

cmin = -max(max(y_mean_zscored(:,:)));
cmax = max(max(y_mean_zscored(:,:)));

cmap = flip(brewermap([],'RdYlBu'));

figure;
s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% Labels and title
xlabel('UMAP1');
ylabel('UMAP2');
title('Frequency');
hold on;

for ii = 1:size(y_mean_zscored,1)
    scatter(reduction(mask==ii,1),reduction(mask==ii,2),10,y_mean_zscored(ii,1:1:length(find(mask==ii))),'filled','MarkerFaceAlpha',1);
    colormap(cmap);
    caxis([cmin cmax])
    caxis([-5 5])
    % caxis([min(min(y_zscore(sst_cell,:))) max(max(y_zscore(sst_cell,:)))])
    % caxis([vmin vmax])
    % caxis([-1 4])
    colorbar;
    hold on;
    % pause;
end

if saveFig
    cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP');
    saveas(s,'vip_spikes_ripples_UMAP.fig')
end

%% Compute structure index

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

% Original emebdding

freq_SI = SI.compute_structure_index(np.array(manifold), np.array(freq), int8(3), params = params );
SI_freq_original = double(freq_SI{1});
SI_freq_shuffled_original = double(freq_SI{4});
SI_freq_99prctile_original = prctile(SI_freq_shuffled_original,99);
SI_freq_95prctile_original = prctile(SI_freq_shuffled_original,95);

amp_SI = SI.compute_structure_index(np.array(manifold), np.array(amp), int8(3), params = params );
SI_amp_original = double(amp_SI{1});
SI_amp_shuffled_original = double(amp_SI{4});
SI_amp_99prctile_original = prctile(SI_amp_shuffled_original,99);
SI_amp_95prctile_original = prctile(SI_amp_shuffled_original,95);

for ii = 1:size(y_mean_zscored,1)
    disp(['Computing vip: ', num2str(ii)]);
    vip_SI = SI.compute_structure_index(np.array(manifold(mask==ii,:)), np.array(y_mean_zscored(ii,1:1:length(find(mask==ii)))), int8(3), params = params );
    SI_vip_original(ii) = double(vip_SI{1});
    SI_vip_shuffled_original = double(vip_SI{4});
    SI_vip_99prctile_original(ii) = prctile(SI_vip_shuffled_original,99);
    SI_vip_95prctile_original(ii) = prctile(SI_vip_shuffled_original,95);
end


% Reduced emebdding

freq_SI = SI.compute_structure_index(np.array(reduction), np.array(freq), int8(3), params = params );
SI_freq_reduction = double(freq_SI{1});
SI_freq_shuffled_reduction = double(freq_SI{4});
SI_freq_99prctile_reduction = prctile(SI_freq_shuffled_reduction,99);
SI_freq_95prctile_reduction = prctile(SI_freq_shuffled_reduction,95);

amp_SI = SI.compute_structure_index(np.array(reduction), np.array(amp), int8(3), params = params );
SI_amp_reduction = double(amp_SI{1});
SI_amp_shuffled_reduction = double(amp_SI{4});
SI_amp_99prctile_reduction = prctile(SI_amp_shuffled_reduction,99);
SI_amp_95prctile_reduction = prctile(SI_amp_shuffled_reduction,95);

for ii = 1:size(y_mean_zscored,1)
    disp(['Computing vip: ', num2str(ii)]);
    vip_SI = SI.compute_structure_index(np.array(reduction(mask==ii,:)), np.array(y_mean_zscored(ii,1:1:length(find(mask==ii)))), int8(3), params = params );
    SI_vip_reduction(ii) = double(vip_SI{1});
    SI_vip_shuffled_reduction = double(vip_SI{4});
    SI_vip_99prctile_reduction(ii) = prctile(SI_vip_shuffled_reduction,99);
    SI_vip_95prctile_reduction(ii) = prctile(SI_vip_shuffled_reduction,95);
end



% Reduction embedding

VIP.SI_original.amp.value = SI_amp_original;
VIP.SI_original.amp.prctile99 = SI_amp_99prctile_original;
VIP.SI_original.amp.prctile95 = SI_amp_95prctile_original;

VIP.SI_original.freq.value = SI_freq_original;
VIP.SI_original.freq.prctile99 = SI_freq_99prctile_original;
VIP.SI_original.freq.prctile95 = SI_freq_95prctile_original;

VIP.SI_original.spikes.value = SI_vip_original;
VIP.SI_original.spikes.prctile99 = SI_vip_99prctile_original;
VIP.SI_original.spikes.prctile95 = SI_vip_95prctile_original;

VIP.SI_reduction.amp.value = SI_amp_reduction;
VIP.SI_reduction.amp.prctile99 = SI_amp_99prctile_reduction;
VIP.SI_reduction.amp.prctile95 = SI_amp_95prctile_reduction;

VIP.SI_reduction.freq.value = SI_freq_reduction;
VIP.SI_reduction.freq.prctile99 = SI_freq_99prctile_reduction;
VIP.SI_reduction.freq.prctile95 = SI_freq_95prctile_reduction;

VIP.SI_reduction.spikes.value = SI_vip_reduction;
VIP.SI_reduction.spikes.prctile99 = SI_vip_99prctile_reduction;
VIP.SI_reduction.spikes.prctile95 = SI_vip_95prctile_reduction;

cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP');
save(['VIP_structure_information_3sessions.mat'],'VIP');
close all;
