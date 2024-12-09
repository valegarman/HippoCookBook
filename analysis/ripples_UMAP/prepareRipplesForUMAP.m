function [ripples,spikes_ripple_UMAP] = prepareRipplesForUMAP(ripples,varargin)


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'bandpass',[80 200]);
addParameter(p,'durations',[-0.05 0.05]);
addParameter(p,'duration_2',[-0.025 0.025]);

parse(p,varargin{:})

basepath = p.Results.basepath;
rippleChannel = p.Results.rippleChannel;
bandpass = p.Results.bandpass;
durations = p.Results.durations;
durations_2 = p.Results.duration_2;

% Care with this line, just a trial
durations = durations_2;

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
% [stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);
% 
% [stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
% [stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
% [stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);

% fastRipple index and entropy (as computed in Foffani et al, 2007; Valero et al, 2017)
% [r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',[-0.5 0.5]);
% matDouble = SyncMap(r,i,'durations',[-0.5 0.5],'smooth',0,'nbins',1250);

% disp('Computing ripple features... ');
% ripSpec = ripSpectrogram(double(matDouble'), samplingRate,0);

% collect results
% data.spectralEntropy = ripSpec.entropyDada;
% data.fastRippleIndex = ripSpec.frippindex;
% data.multiTapperFreq = ripSpec.freqData;
% 
% maps.multitaperSpecs = ripSpec;

stats = [];
rippleStats.stats = stats;
rippleStats.data = data;
rippleStats.maps = maps;

ripples.stats = stats;
ripples.data = data;
ripples.maps = maps;

ripples.rippleStats = rippleStats;

ripples.maps_filtered = maps.ripples_filtered;

% Compute pairwise distances
D = pdist2(ripples.maps_filtered, ripples.maps_filtered);

% Set diagonal elements to NaN
D(1:size(D, 1)+1:end) = NaN;

% Compute nearest neighbor distance threshold (5th percentile)
nn_dist = sum(D < prctile(D(:), 5), 2);

% Find noise points based on nearest neighbor distances (20th percentile)
noiseIdx = nn_dist < prctile(nn_dist, 20);

ripples.maps_filtered(find(noiseIdx),:) = [];
ripples.peaks(find(noiseIdx)) = [];
ripples.timestamps(find(noiseIdx),:) = [];

ripples.data.peakFrequency(find(noiseIdx),:) = [];
ripples.data.peakAmplitude(find(noiseIdx),:) = [];

save(['ripple_UMAP.mat'],'ripples')

ripples_UMAP = ripples.maps_filtered;

ripples_UMAP_v2.maps_filtered = ripples_UMAP;
ripples_UMAP_v2.peakFrequency = ripples.data.peakFrequency;
ripples_UMAP_v2.peakAmplitude = ripples.data.peakAmplitude;

save(['UMAP_ripples_v2.mat'],'ripples_UMAP_v2');

% Compute firing rate during ripple window (+- 0.025)

spikes_ripple_UMAP = nan(length(spikes.UID),2142);

for ii = 1:length(spikes.UID)
    parfor jj = 1:length(ripples.peaks)

        spikes_ripple_UMAP_v2(ii,jj) = length(find(InIntervals(spikes.times{ii},[ripples.peaks(jj)+durations(1) ripples.peaks(jj)+durations(2)])));
        
    end
end

save('spikes_ripples_UMAP_v2.mat','spikes_ripple_UMAP_v2');

end