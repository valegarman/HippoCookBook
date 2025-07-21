
close all;
clear all;

%% Load the LFP file 
session = loadSession;
Channels = session.extracellular.electrodeGroups.channels{3};
lfp = getLFP(Channels,'noPrompts', true);   % lfp datas
lfp_data = double(lfp(1).data');            % Channels x Time sampling


%% pre-processing the signal, filtering and sampling

Fs = lfp(1).samplingRate;        % sampling rate 
nChannels = length(Channels);    % n channels 

% design the filter FIR bandpass 
thetaFilt = designfilt('bandpassfir','FilterOrder', 100,'CutoffFrequency1', 5,'CutoffFrequency2', 12,'SampleRate', Fs);
slowGammaFilt = designfilt('bandpassfir', 'FilterOrder', 100,'CutoffFrequency1', 30,'CutoffFrequency2', 60,'SampleRate', Fs);
fastGammaFilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1', 60,'CutoffFrequency2', 90,'SampleRate', Fs);
rippleFilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1', 100,'CutoffFrequency2', 250,'SampleRate', Fs);

% design the matrix for the data filtered in each band 
theta_data = zeros(nChannels, size(lfp_data,2));
slowgamma_data = zeros(nChannels, size(lfp_data,2));
fastgamma_data = zeros(nChannels, size(lfp_data,2));
ripple_data = zeros(nChannels, size(lfp_data,2));

% filtering and perform the hilbert transformation to get the envolope for each band 
for ch = 1:nChannels
    raw = lfp_data(ch,:); % signal for each channel 
    theta_data(ch,:)      = abs(hilbert(filtfilt(thetaFilt, raw)));
    slowgamma_data(ch,:)  = abs(hilbert(filtfilt(slowGammaFilt, raw)));
    fastgamma_data(ch,:)  = abs(hilbert(filtfilt(fastGammaFilt, raw)));
    ripple_data(ch,:)     = abs(hilbert(filtfilt(rippleFilt, raw)));
end

% Generate the feature vector [time × (channels × 4)]
features = [theta_data', slowgamma_data', fastgamma_data', ripple_data'];

features = [theta_data', ripple_data'];

%% temporal segmentation to decrease the computational cost 

bin_size = 0.5;                                   % 0.5 seconds
bin_samples = round(bin_size * Fs);               % samples per bin      
nBins = floor(size(features,1) / bin_samples);
        
binned_features = zeros(nBins, size(features,2)); % [nBins x (nChannels x 4)]

for i = 1:nBins
    idx_start = (i-1)*bin_samples + 1;
    idx_end = idx_start + bin_samples - 1;
    binned_features(i,:) = mean(features(idx_start:idx_end, :), 1); % media per bin
end

binned_z = zscore(binned_features);   % Normalizing 
binned_z = binned_z(1:2:end, :);      % Downsampling - taking only 1 each 2 points

%% The paper used a non linear dimensionality reduction technic PHATE to obtain a low dimension embedding of the filtered LFP 
% PHATE (Potential of Heat-diffusion for Affinity-based Transition Embedding) 

%% 1. t-SNE reduction
perplexity = 50;
n_dims = 2;

% Applica t-SNE
tsne_embedding = tsne(binned_z, 'NumDimensions', n_dims, 'Perplexity', perplexity);

% Visualizza
figure;
scatter(tsne_embedding(:,1), tsne_embedding(:,2), 15, linspace(1, size(tsne_embedding,1), size(tsne_embedding,1)), 'filled');
colormap(parula);
colorbar;
title('t-SNE embedding of LFP features');
xlabel('t-SNE 1');
ylabel('t-SNE 2');

%% 2. UMAP
% Applica UMAP
[umap_embedding, umap_params] = run_umap(binned_z, 'n_components', 2);

% Visualizza
figure;
scatter(umap_embedding(:,1), umap_embedding(:,2), 15, linspace(1, size(umap_embedding,1), size(umap_embedding,1)), 'filled');
colormap(parula);
colorbar;
title('UMAP embedding of LFP features');
xlabel('UMAP 1');
ylabel('UMAP 2');

%% 3. PCA
% PCA
[coeff, score, latent] = pca(binned_z);

% Prendi le prime 3 componenti per visualizzare il manifold
manifold3D = score(:,1:3); % [nBins x 3]

% Visualizzazione del manifold
figure;
scatter3(manifold3D(:,1), manifold3D(:,2), manifold3D(:,3), 10, 'filled');
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
title('Manifold dello stato di rete (PCA)');

% t-SNE
Y = tsne(binned_z, 'NumDimensions', 2, 'Perplexity', 30);
figure;
scatter(Y(:,1), Y(:,2), 10, 'filled');
title('Manifold t-SNE');
addpath('drtoolbox/techniques/');

%% 4. ISOMAP with geodetic distance and k-nearest neighbors

% Dati: binned_z [nBins x nFeatures], normalizzato (zscore)
% Calcola distanza euclidea
D = pdist2(binned_z, binned_z);

% k-nearest neighbors (Isomap richiede un grafo)
k = 10;
n = size(D,1);
for i = 1:n
    [~, idx] = sort(D(i,:));
    D(i, idx(k+2:end)) = Inf;  % tiene solo i k vicini
end


% Rendi la matrice simmetrica (grafo non orientato)
D = min(D, D');

% Calcola distanze geodetiche (Floyd–Warshall)
G = graph(D);
geoD = distances(G);  % distanza lungo il grafo

% MDS (Multidimensional Scaling) per embedding
Y = mdscale(geoD, 2);  % proietta in 2D

figure;
scatter(Y(:,1), Y(:,2), 15, linspace(1,n,n), 'filled');
title('Isomap embedding of LFP features');
xlabel('Dim 1');
ylabel('Dim 2');
colormap(parula);
colorbar;

% Applica Isomap
[Y, R] = compute_mapping(binned_z, 'Isomap', 2, k);  % 2D embedding, k vicini

% Visualizza
figure;
scatter(Y(:,1), Y(:,2), 15, linspace(1, size(Y,1), size(Y,1)), 'filled');
title('Isomap embedding');

% design the filter FIR bandpass 
Wp_theta = [5 12] / (Fs/2);  
[b_theta, a_theta] = cheby1(4, 0.5, Wp_theta, 'bandpass');

Wp_sg = [30 60] / (Fs/2);
[b_sg, a_sg] = cheby1(4, 0.5, Wp_sg, 'bandpass');

Wp_fg = [60 90] / (Fs/2);
[b_fg, a_fg] = cheby1(4, 0.5, Wp_fg, 'bandpass');

Wp_rip = [100 250] / (Fs/2);
[b_rip, a_rip] = cheby1(4, 0.5, Wp_rip, 'bandpass');

% design the matrix for the data filtered in each band 
theta_data = zeros(nChannels, size(lfp_data,2));
slowgamma_data = zeros(nChannels, size(lfp_data,2));
fastgamma_data = zeros(nChannels, size(lfp_data,2));
ripple_data = zeros(nChannels, size(lfp_data,2));

% filtering and perform the hilbert transformation to get the envolope for each band 
for ch = 1:nChannels
    raw = lfp_data(ch,:); % signal for each channel 
    theta_data(ch,:) = abs(hilbert(filtfilt(b_theta, a_theta, raw)));
    slowgamma_data(ch,:) = abs(hilbert(filtfilt(b_sg, a_sg, raw)));
    fastgamma_data(ch,:) = abs(hilbert(filtfilt(b_fg, a_fg, raw)));
    ripple_data(ch,:)    = abs(hilbert(filtfilt(b_rip, a_rip, raw)));
end


% For SWR event detection, the LFPs were initially referenced against 952 a channel where CA1 ripples were not observed.
ripple_band = [130 200]; % ripple band
powerProfile = powerSpectrumProfile(ripple_band, 'showfig', false,'saveMat',false);

[~, idxMax] = max(powerProfile.mean);
maxRippleChannel = powerProfile.channels(idxMax);
[ripples] = rippleMasterDetector('rippleChannel', maxRippleChannel,'SaveMat',false);

signal_interest = lfp_data(maxRippleChannel, :);
signal_ref = lfp_data(reference_channel, :);  % reference_channel va scelto come canale senza ripple

lfp_diff = signal_interest - signal_ref;


% Supponiamo che 'powerProfile' abbia anche la potenza ripple per tutti i canali
ripplePower = powerProfile.mean; % potenza media in banda ripple

% Trova canale con potenza ripple minima (escludendo maxRippleChannel)
candidateChannels = setdiff(1:nChannels, maxRippleChannel);
[~, idxRef] = min(ripplePower(candidateChannels));
reference_channel = candidateChannels(idxRef);















% Generate the feature vector [time × (channels × 4)]
features = [theta_data', slowgamma_data', fastgamma_data', ripple_data'];

features = [theta_data', ripple_data'];

%% temporal segmentation to decrease the computational cost 

bin_size = 0.5;                                   % 0.5 seconds
bin_samples = round(bin_size * Fs);               % samples per bin      
nBins = floor(size(features,1) / bin_samples);
        
binned_features = zeros(nBins, size(features,2)); % [nBins x (nChannels x 4)]

for i = 1:nBins
    idx_start = (i-1)*bin_samples + 1;
    idx_end = idx_start + bin_samples - 1;
    binned_features(i,:) = mean(features(idx_start:idx_end, :), 1); % media per bin
end

binned_z = zscore(binned_features);   % Normalizing 
binned_z = binned_z(1:2:end, :);      % Downsampling - taking only 1 each 2 points

%% The paper used a non linear dimensionality reduction technic PHATE to obtain a low dimension embedding of the filtered LFP 
% PHATE (Potential of Heat-diffusion for Affinity-based Transition Embedding) 

%% 1. t-SNE reduction
perplexity = 50;
n_dims = 2;

% Applica t-SNE
tsne_embedding = tsne(binned_z, 'NumDimensions', n_dims, 'Perplexity', perplexity);

% Visualizza
figure;
scatter(tsne_embedding(:,1), tsne_embedding(:,2), 15, linspace(1, size(tsne_embedding,1), size(tsne_embedding,1)), 'filled');
colormap(parula);
colorbar;
title('t-SNE embedding of LFP features');
xlabel('t-SNE 1');
ylabel('t-SNE 2');

%% 2. UMAP
% Applica UMAP
[umap_embedding, umap_params] = run_umap(binned_z, 'n_components', 2);

% Visualizza
figure;
scatter(umap_embedding(:,1), umap_embedding(:,2), 15, linspace(1, size(umap_embedding,1), size(umap_embedding,1)), 'filled');
colormap(parula);
colorbar;
title('UMAP embedding of LFP features');
xlabel('UMAP 1');
ylabel('UMAP 2');

%% 3. PCA
% PCA
[coeff, score, latent] = pca(binned_z);

% Prendi le prime 3 componenti per visualizzare il manifold
manifold3D = score(:,1:3); % [nBins x 3]

% Visualizzazione del manifold
figure;
scatter3(manifold3D(:,1), manifold3D(:,2), manifold3D(:,3), 10, 'filled');
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
title('Manifold dello stato di rete (PCA)');

% t-SNE
Y = tsne(binned_z, 'NumDimensions', 2, 'Perplexity', 30);
figure;
scatter(Y(:,1), Y(:,2), 10, 'filled');
title('Manifold t-SNE');
addpath('drtoolbox/techniques/');

%% 4. ISOMAP with geodetic distance and k-nearest neighbors

% Dati: binned_z [nBins x nFeatures], normalizzato (zscore)
% Calcola distanza euclidea
D = pdist2(binned_z, binned_z);

% k-nearest neighbors (Isomap richiede un grafo)
k = 10;
n = size(D,1);
for i = 1:n
    [~, idx] = sort(D(i,:));
    D(i, idx(k+2:end)) = Inf;  % tiene solo i k vicini
end


% Rendi la matrice simmetrica (grafo non orientato)
D = min(D, D');

% Calcola distanze geodetiche (Floyd–Warshall)
G = graph(D);
geoD = distances(G);  % distanza lungo il grafo

% MDS (Multidimensional Scaling) per embedding
Y = mdscale(geoD, 2);  % proietta in 2D

figure;
scatter(Y(:,1), Y(:,2), 15, linspace(1,n,n), 'filled');
title('Isomap embedding of LFP features');
xlabel('Dim 1');
ylabel('Dim 2');
colormap(parula);
colorbar;

% Applica Isomap
[Y, R] = compute_mapping(binned_z, 'Isomap', 2, k);  % 2D embedding, k vicini

% Visualizza
figure;
scatter(Y(:,1), Y(:,2), 15, linspace(1, size(Y,1), size(Y,1)), 'filled');
title('Isomap embedding');
