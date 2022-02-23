function [stats] = computeRippleStats(varargin)

%RippleStats - Compute descriptive stats for ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples,<options>)
%
%    filtered       ripple-band filtered samples (one channel)
%    timestamps	    timestamps (in seconds) to match filtered vector
%    ripples        ripple timing information (STRUCT VERSION) (obtained using <a href="matlab:help bz_FindRipples">bz_FindRipples</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bandpass'    Default, [80 200]
%     'ripples'     By default tries to load ripple.event.mat structure or
%                   runs findRipples. Otherwise, ripple timing information structure(see
%                       findRipples)
%     'durations'   durations before and after ripple peak (in s)
%                   (default = [-0.075 0.075])
%     'basepath'    Default, pwd
%     'saveSummary' Default, true
%     'saveMat'     Save or updates ripples.event.mat structure with status info.
%                       Default, true
%     'plotOpt'     Default, true
%     'rippleChannel'
%                   By default loads getHippocampalLayers, and uses
%                   hippocampalLayer.bestShankLayers.pyramidal (1-index)
%    =========================================================================
%
%  OUTPUT
%
%    maps.ripples               instantaneous amplitude (one ripple per row)
%    maps.frequency             instantaneous frequency
%    maps.phase                 instantaneous phase
%    maps.amplitude             enveloppe amplitude
%    data.peakFrequency         frequency at peak
%    data.peakAmplitude         amplitude at peak
%    data.duration              durations
%    stats.durationAmplitude    correlation between duration and amplitude (rho, p)
%    stats.durationFrequency    correlation between duration and frequency (rho, p)
%    stats.amplitudeFrequency   correlation between amplitude and frequency (rho, p)
%    stats.acg.data             autocorrelogram data
%    stats.acg.t                autocorrelogram time bins
%    and more...

% edited by David Tingley to fit buzcode formatting standards, 2017
% edited by Manuel Valero 2022 (fix inputs, include FRindex and entropy, as in Foffani et al, 2007; Valero et al, 2017)
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'bandpass',[80 200],@isnumeric);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'durations',[-0.075 0.075], @isnumeric);
addParameter(p,'ripples',[], @isstruct);
addParameter(p,'rippleChannel',[], @isstruct);

parse(p,varargin{:})
basepath = p.Results.basepath;
bandpass = p.Results.bandpass;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
durations = p.Results.durations;
ripples = p.Results.ripples;
rippleChannel = p.Results.rippleChannel;

% Dealing with inputs
prevBasepath = pwd;
cd(basepath);

% session = sessionTemplate(basepath,'showGUI',false);
if isempty(rippleChannel)
    hippocampalLayers = getHippocampalLayers;
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

if isempty(ripples)
    ripples = findRipples(rippleChannel,'thresholds',[1.5 3.5],'passband',bandpass,...
        'EMGThresh',1,'durations',[20 250], 'saveMat', saveMat);
end

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

% idx = ceil((ripples.peaks-ripples.timestamps(:,1))*ripples.detectorinfo.detectionparms.frequency);

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

maps.timestamps = linspace(durations(1),durations(2),nBins)

data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = abs(diff(ripples.timestamps'))';

% stats
[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);
[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);

% fastRipple index and entropy (as computed in Foffani et al, 2007; Valero et al, 2017)
keyboard;
[r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',[-0.5 0.5]);
matDouble = SyncMap(r,i,'durations',[-0.5 0.5],'smooth',0,'nbins',1250);

ripSpec = ripSpectrogram(double(matDouble), samplingRate,0)
figure;


resolHz=10; %spectral resolution in Hz
nfft=floor(samplingRate/resolHz);
nw=2; %parametro del multitaper ("time-bandwidth product")
noverlap=nfft/1.024; %window overlap default:nfft/2
aa=(100/resolHz)+1; bb=(600/resolHz)+1; cc=(400/resolHz)+1; %aa=11; bb=61; cc=41;
for ii = 1:size(maps.ripples_raw,1) 
    for jj = 1:size(maps.ripples_raw,2)
        [S,F]= pmtm(maps.ripples_raw(ii,:),nw,nfft,samplingRate);
    end
end
T=(1:floor(length(X)/delta)*delta)*(1/Fs);

[TimeFreq(:,:,i),Ftf,T] = SpectrogramMultiTaper(maps.ripples_raw,noverlap,nfft,samplingRate,nw);


keyboard;

rippleStats.stats = stats;
rippleStats.data = data;
rippleStats.maps = maps;


end
