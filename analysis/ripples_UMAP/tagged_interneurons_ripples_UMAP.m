
clearvars -except projectResults projectSessionResults taggedCells 
% cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data')
% file = dir('chapter_0_data19-Sep-2024.mat');
% load(file.name);

% General configuration
addpath(genpath('C:\Users\pabad\Downloads\umapFileExchange (4.4)\umap'));
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

verbose = true;
cmap = flip(brewermap([],'RdYlBu'));

% Dimensionality reduction parameters
n_neighbors_UMAP = 15;
min_dist = 0.1;
D = 4; % (Intrinsic dimension based on Liset's paper)
metric = 'euclidean';

genotype = {'pv','sst','id2','vip','camk2'};

for ii = 1:length(genotype)

    units.(genotype{ii}) = taggedCells.hippo.(genotype{ii});
    sessions.(genotype{ii}) = unique(projectResults.session(taggedCells.hippo.(genotype{ii})))

    for jj = 1:length(sessions.(genotype{ii}))

        cell = find(taggedCells.hippo.(genotype{ii})(ismember(projectResults.session,lower(sessions.(genotype{ii}){jj}))))

        animal = strsplit(sessions.(genotype{ii}){jj},'_');
        animal = animal{1};
        % Loading ripples_UMAP and spikes_ripples_UMAP
        try
            cd(['Z:\',animal,'\',sessions.(genotype{ii}){jj}])
        catch
            cd(['Y:\',animal,'\',sessions.(genotype{ii}){jj}])
        end

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
        SI_amp_original.(genotype{ii}).value(jj) = double(amp_SI{1});
        SI_amp_shuffled_original = double(amp_SI{4});
        SI_amp_original.(genotype{ii}).prctile99(jj) = prctile(SI_amp_shuffled_original,99);
        SI_amp_original.(genotype{ii}).prctile95(jj) = prctile(SI_amp_shuffled_original,95);

        amp_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(amp), int8(3), params = params );
        SI_amp_reduction.(genotype{ii}).value(jj) = double(amp_reduction_SI{1});
        SI_amp_reduction_shuffled = double(amp_reduction_SI{4});
        SI_amp_reduction.(genotype{ii}).prctile99(jj) = prctile(SI_amp_reduction_shuffled,99);
        SI_amp_reduction.(genotype{ii}).prctile95(jj) = prctile(SI_amp_reduction_shuffled,95);

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
        SI_freq_original.(genotype{ii}).value(jj) = double(freq_SI{1});
        SI_freq_shuffled_original = double(freq_SI{4});
        SI_freq_original.(genotype{ii}).prctile99(jj) = prctile(SI_freq_shuffled_original,99);
        SI_freq_original.(genotype{ii}).prctile95(jj) = prctile(SI_freq_shuffled_original,95);

        freq_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(freq), int8(3), params = params );
        SI_freq_reduction.(genotype{ii}).value(jj) = double(freq_reduction_SI{1});
        SI_freq_reduction_shuffled = double(freq_reduction_SI{4});
        SI_freq_reduction.(genotype{ii}).prctile99(jj) = prctile(SI_freq_reduction_shuffled,99);
        SI_freq_reduction.(genotype{ii}).prctile95(jj) = prctile(SI_freq_reduction_shuffled,95);

        % Now spikes projected

        for kk = 1:length(cell)

            y_zscore = zscore(spikes_UMAP_v2(cell(kk),:));

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
            SI_spikes_original.(genotype{ii}).(sessions.(genotype{ii}){jj}).value(kk) = double(spikes_SI{1}); 
            SI_spikes_shuffled_original = double(spikes_SI{4});
            SI_spikes_original.(genotype{ii}).(sessions.(genotype{ii}){jj}).prctile99(kk) = prctile(SI_spikes_shuffled_original,99);
            SI_spikes_original.(genotype{ii}).(sessions.(genotype{ii}){jj}).prctile95(kk) = prctile(SI_spikes_shuffled_original,95);
    
            spikes_reduction_SI = SI.compute_structure_index(np.array(reduction), np.array(y_zscore), int8(3), params = params );
            SI_spikes_reduction.(genotype{ii}){jj}.value(kk) = double(spikes_reduction_SI{1});
            SI_spikes_reduction_shuffled = double(spikes_reduction_SI{4});
            SI_spikes_reduction.(genotype{ii}){jj}.prctile99(kk) = prctile(SI_spikes_reduction_shuffled,99);
            SI_spikes_reduction.(genotype{ii}){jj}.prctile95(kk) = prctile(SI_spikes_reduction_shuffled,95);

            % disp(SI_spikes_original.(genotype{ii}).(sessions.(genotype{ii}){jj}).value(kk));
            % disp(SI_spikes_reduction.(genotype{ii}){jj}.value(kk));

        end
    end


end

% Save the data

SI = [];

cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data\UMAP ripples');

save('tagged_SI_freq_original.mat','SI_freq_original');
save('tagged_SI_freq_reduction.mat','SI_freq_reduction');

save('tagged_SI_amp_original.mat','SI_amp_original');
save('tagged_SI_amp_reduction.mat','SI_amp_reduction');


save('tagged_SI_spikes_original.mat','SI_spikes_original');
save('tagged_SI_spikes_reduction.mat','SI_spikes_reduction');


% Plotting
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

flds = fields(SI_spikes_original);

SI = [];


for ii = 1:length(flds)
    
    SI.(flds{ii}).value = [];
    SI.(flds{ii}).p99 = [];
    SI.(flds{ii}).p95 = [];

    sessions = fields(SI_spikes_original.(flds{ii}));

    for jj = 1:length(sessions)

        for kk = 1:length(SI_spikes_original.(flds{ii}).(sessions{jj}).value)
            SI.(flds{ii}).value = [SI.(flds{ii}).value SI_spikes_original.(flds{ii}).(sessions{jj}).value(kk)];
            SI.(flds{ii}).p99 = [SI.(flds{ii}).p99 SI_spikes_original.(flds{ii}).(sessions{jj}).prctile99(kk)];
            SI.(flds{ii}).p95 = [SI.(flds{ii}).p95 SI_spikes_original.(flds{ii}).(sessions{jj}).prctile95(kk)];
        end
    end
end

figure;

[gs] = groupStats({SI.camk2.value, SI.pv.value, SI.sst.value,...
    SI.vip.value, SI.id2.value},[], ...
    'color',[colors.camk2; colors.pv; colors.sst; colors.vip; colors.id2],'plotData',true,'sigStar',true);

% MAX PV
% fPV4_sess2 (neuron # 18)

cd('Z:\fPv4\fPv4_210309_sess2')

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

y_zscore = zscore(spikes_UMAP_v2(18,:));

if verbose
    figure;
    s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Spikes');
    hold on;
    scatter(reduction(:,1),reduction(:,2),10,y_zscore,'filled','MarkerFaceAlpha',1);
    caxis([-2 2])
    colormap(cmap);
    colorbar;

    cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP')

    saveas(s,['bestPV_rippleUMAP.fig'])
end

%%
% MAX ID2

[a,b] = max(SI.id2.value);
% fnkx9_sess9 (neuron # 8)

cd('Z:\fNkx9\fNkx9_200827_sess9')

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

y_zscore = zscore(spikes_UMAP_v2(8,:));

if verbose
    figure;
    s = scatter(reduction(:, 1), reduction(:, 2), 10, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Spikes');
    hold on;
    scatter(reduction(:,1),reduction(:,2),10,y_zscore,'filled','MarkerFaceAlpha',1);
    caxis([-2 2])
    colormap(cmap);
    colorbar;

    cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\figures\Figures ripples UMAP')

    saveas(s,['bestId2_rippleUMAP.fig'])
end

