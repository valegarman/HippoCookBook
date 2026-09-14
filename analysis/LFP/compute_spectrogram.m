
function [spectogram] = compute_spectrogram(varargin)
% Detect theta/delta periods
% 
% INPUTS
% <optional>
% 'basepath'            Default pwd
% 'lfp'                 buzcode-formatted lfp structure (use bz_GetLFP)
%                           needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%                           If empty or no exist, look for lfp in basePath folder
% 'saveSummary'         Default true
% 'saveMat'             Detault true
% 'force'               Default false
% 'bandpass'            Default [6 12]
% 'powerThreshold'      Default 1 SD
% 'channel'             Numeric [ex, 5]; by default calls
%                           getHippocampalLayers and uses oriens.
% 'updateSleepStates'   Default true
% 'useCSD'              Default, true.
% 'discardRipples'      Discard ripples from nonTheta, default true.
% 
% OUTPUT
% thetaEpochs           states structure with theta epochs intervals
%
% Manu Valero 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'bandpass',[6 12], @isnumeric);
addParameter(p,'channel',[],@isnumeric);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'useCSD',false,@islogical);
addParameter(p,'uselog10Power',true,@islogical);
addParameter(p,'fpass',[2 120],@isnumeric);
addParameter(p,'intervals',[],@isnumerictype);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
bandpass = p.Results.bandpass;
channel = p.Results.channel;
plotting = p.Results.plotting;
useCSD = p.Results.useCSD;
uselog10Power = p.Results.uselog10Power;
fpass = p.Results.fpass;
intervals = p.Results.intervals;

% Deal with inputs
prevBasepath = pwd;
cd(basepath);

% targetFile = dir('*.spectrogram.channelinfo.mat');
% if ~isempty(targetFile) && ~force
%     disp('Theta epochs already detected! Loading file.');
%     load(targetFile.name);
%     return
% end

if isempty(intervals)
    intervals = [0 Inf];
end

if isempty(lfp) && ~useCSD
    lfpT = getLFP(channel,'noPrompts',true, 'intervals', intervals);
elseif useCSD
    disp('Computing CSD...');
    lfpT = computeCSD(lfp,'channels',channel);
else
    warning('CSD estimation not possible. Using LFP...');
end

samplingRate = lfpT.samplingRate;
[wave,f,t,~,wphases,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,bandpass(1),bandpass(2),8,0);
[~,mIdx]=max(wave); % get index max power for each timepiont
pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))]; %converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
lfpphase=wphases(pIdx); %get phase of max amplitude wave at each timepoint
lfpphase = mod(lfpphase,2*pi); %covert to 0-2pi rather than -pi:pi
power = rms(abs(wave))';

if uselog10Power
    power = log10(power);
end

params.Fs = lfpT.samplingRate; params.fpass = fpass; params.tapers = [3 5]; params.pad = 1;
[S,t,f] = mtspecgramc_fast(single(lfpT.data),[2 1],params); S(S==0) = NaN;
S = log10(S); % in Db
%S_det= bsxfun(@minus,S,polyval(polyfit(f,nanmean(S,1),2),f)); % detrending
S_det= detrend(S',2)';

spectogram.lfpphase = lfpphase;
spectogram.samplingRate = samplingRate;
spectogram.power = power;
spectogram.timestamps_power = lfpT.timestamps;
spectogram.timestamps_S = t;
spectogram.fpass = fpass;
spectogram.tapers = params.tapers;
spectogram.frequencies = f;
spectogram.S = S;
spectogram.bandpass = bandpass;
spectogram.S_det = S_det;
% spectogram.intervals = intervals;
spectogram.channel = channel;
spectogram.params.bandpass = bandpass;
spectogram.params.uselog10Power = uselog10Power;
spectogram.params.useCSD = useCSD;

if saveMat
    disp('Saving...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.spectrogram.channelifo_' num2str(channel) 'ch.mat'],'spectogram');
end

if plotting

    figure;
    subplot(3,3,[1 2 4 5])
    imagesc(t/60,f,S_det',[-1.5 1.5]);
    ylim([1.5 50]);
    set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('');
    title(['Channel: ' num2str(channel)],'FontWeight','normal');
    ylim([2.5 50]);

    subplot(3,3,[7 8])
    plot(spectogram.timestamps_power/60, smooth(power, samplingRate*10), 'k');
    xlabel('Time [min]');
    set(gca,'TickDir','out');
    ylabel(['Power (dB,' num2str(bandpass(1)) '-' num2str(bandpass(2)) ')']);
    colormap jet
    
    subplot(1,3,[3])
    plotFill(f,S,'color', [.3 .3 .3],'lineStyle', '-', 'xscale', 'linear'); xlim([1 120]);
    ylabel('Power [dB]'); xlabel('Freq [Hz]');   
    
    if saveSummary
        mkdir('SummaryFigures'); % create folder
        saveas(gcf,'SummaryFigures\spectrogram.png');
    end
end


cd(prevBasepath);
end