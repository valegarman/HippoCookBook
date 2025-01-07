function [UMAP,spikes_ripple_UMAP] = getRipplesUMAP(ripples,varargin)


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'bandpass',[80 200]);
addParameter(p,'durations',[-0.025 0.025]);

parse(p,varargin{:})

basepath = p.Results.basepath;
rippleChannel = p.Results.rippleChannel;
bandpass = p.Results.bandpass;
durations = p.Results.durations;

% Load ripple.mat
% file = dir('*.ripplesAll.event.mat');
% load(file.name);

rippleChannel = ripples.detectorinfo.detectionchannel;

% Load spikes
spikes = loadSpikes();

% some fixed parameters
corrBinSize = 0.01;

% get filtered signal
lfp = getLFP(rippleChannel);
filteredLFP = bz_Filter(lfp,'channels',rippleChannel,'filter','butter','passband',bandpass,'order',3);
samplingRate = filteredLFP.samplingRate;
timestamps = filteredLFP.timestamps;
filtered = filteredLFP.data;
unfiltered = lfp.data;

nBins = floor(samplingRate*diff(durations)/2)*2+1; % must be odd
nHalfCenterBins = 3;
centerBin = ceil(nBins/2);
centerBins = centerBin-nHalfCenterBins:centerBin+nHalfCenterBins;

% Get instantaneous phase and amplitude
phase = filteredLFP.phase;
amplitude = abs(filteredLFP.hilb);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),timestamps,'smooth',0);
frequency = frequency/(2*pi);

% Compute ripple map
[r,i] = Sync([timestamps filtered],ripples.peaks,'durations',durations);
maps.ripples_filtered = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);
[r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',durations);
maps.ripples_raw = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute frequency Map
[f,i] = Sync([timestamps frequency],ripples.peaks,'durations',durations);
maps.frequency = SyncMap(f,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute phase map
[p,i] = Sync([timestamps phase],ripples.peaks,'durations',durations);
maps.phase = SyncMap(p,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute amplitude map
[a,i] = Sync([timestamps amplitude],ripples.peaks,'durations',durations);
maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

maps.timestamps = linspace(durations(1),durations(2),nBins);

data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = abs(diff(ripples.timestamps'))';

% InterRipple Frequency
data.interRippleFrequency = [NaN; 1./diff(ripples.peaks)];

% stats
[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);

[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);

% fastRipple index and entropy (as computed in Foffani et al, 2007; Valero et al, 2017)
[r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',[-0.5 0.5]);
matDouble = SyncMap(r,i,'durations',[-0.5 0.5],'smooth',0,'nbins',1250);

disp('Computing ripple features... ');
ripSpec = ripSpectrogram(double(matDouble'), samplingRate,0);

% collect results
data.spectralEntropy = ripSpec.entropyDada;
data.fastRippleIndex = ripSpec.frippindex;
data.multiTapperFreq = ripSpec.freqData;

maps.multitaperSpecs = ripSpec;

rippleStats.stats = stats;
rippleStats.data = data;
rippleStats.maps = maps;

ripples.stats = stats;
ripples.data = data;
ripples.maps = maps;

ripples.rippleStats = rippleStats;

ripples.maps_filtered_raw = ripples.maps.ripples_filtered;

save(['ripple_UMAP.mat'],'ripples')

% Compute firing rate during ripple window (+- 0.025)

for ii = 1:length(spikes.UID)
    for jj = 1:length(ripples.peaks)

        spikes_ripple_UMAP(ii,jj) = length(find(InIntervals(spikes.times{ii},[ripples.peaks(jj)-0.02 ripples.peaks(jj)+0.025])));
        
    end
end

save('spikes_ripples_UMAP.mat','spikes_ripple_UMAP');

%% UMAP

manifold = ripples.rippleStats.maps.ripples_filtered;

% Dimensionality reduction

n_neighbors = 15;
min_dist = 0.1;
D = 2;
metric = 'euclidean';

[reduction,umap,clusterIdentifiers,extras] = run_umap(manifold,...
    'min_dist',min_dist,'n_neighbors',n_neighbors,'metric',metric,...
    'n_components',D);

figure;
% ax = subplot(1, 3, 1);
% Scatter plot
s = scatter(reduction(:, 1), reduction(:, 2), 50, ripples.rippleStats.data.peakFrequency, 'filled', 'MarkerEdgeColor', 'k');
% colormap(ax, 'coolwarm');
colormap(gcf, 'pink');
caxis([prctile(ripples.rippleStats.data.peakFrequency, 5), prctile(ripples.rippleStats.data.peakFrequency, 95)]); % Set color limits

% Labels and title
xlabel('UMAP1');
ylabel('UMAP2');
title('Frequency');
% Colorbar
colorbar('eastoutside');

UMAP = [];
UMAP.reduction = reduction;
% UMAP.umap = umap;
UMAP.clusterIdentifiers = clusterIdentifiers;
% UMAP.extras = extras;



clear umap
save('UMAP_for_ripples.mat','UMAP');


figure;
% ax = subplot(1, 3, 1);
% Scatter plot
s = scatter(reduction(:, 1), reduction(:, 2), 50, ripples.rippleStats.data.peakFrequency, 'filled', 'MarkerEdgeColor', 'k');
% colormap(ax, 'coolwarm');
colormap(gcf, 'pink');
caxis([prctile(ripples.rippleStats.data.peakFrequency, 5), prctile(ripples.rippleStats.data.peakFrequency, 95)]); % Set color limits


end