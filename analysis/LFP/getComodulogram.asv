function [ powerComod ] = getComodulogram(lfp,varargin)
%[ powerComod ] = bz_MyComodulogram(lfp,varargin) calculates the
%comodulogram (i.e. power-power correlation) for an lfp file.
%
%
%INPUT
%   LFP         LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps, lfp.refCh
%               (optional). If only one lfp channel is provided, that is unsed
%               for corr computation. If two, pairwise correlation is
%               performed. If one refCh and n data channels are provided
%               the correlation of the refCh against all data channels is
%               performed. 
% <OPTIONAL>
%
%   frange      [min max] frequency
%   nfreqs      number of frequencies 
%   space       spacing between freqs, 'lin' or 'log'
%   fvector     predefined vector of frequencies 
%   specnorm    normalization for spectral power,
%                   options: 'mean','logmean','log' (default: log)
%   numvarbins  number of bins for your external variable
%   varnorm     normalization for the variable,
%                   options: 'percentile', 'none' (default: 'none')
%   type        'wavelet' or 'FFT'      (default: 'wavelet')
%      -if type: 'wavelet'-
%       .ncyc       number of wavelet cycles (recommended: ~5)
%      -if type: 'FFT'-
%       .winsize (s)
%       .overlap
%   doPlot      Default true
%   saveMat     Default true
%
% DLevenstein 2017
% Modified by Antonio FR, 7/18/18
% Parsong modified by Manu Valero, 5/9/2020
%
%
%% Parse the inputs
%Parameters
p = inputParser;
addParameter(p,'frange',[1 128],@isnumeric);
addParameter(p,'nfreqs',100,@isnumeric);
addParameter(p,'ncyc',5,@isnumeric);
addParameter(p,'space','log');
addParameter(p,'samplingRate',[]);
addParameter(p,'showprogress',true,@islogical);
addParameter(p,'saveMat',false);
addParameter(p,'fvector',[]);
addParameter(p,'specnorm','log');
addParameter(p,'type','wavelet');
addParameter(p,'winsize',1,@isnumeric);
addParameter(p,'overlap',0.5,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'force',false,@islogical)

parse(p,varargin{:});
specparms.frange = p.Results.frange;
specparms.nfreqs = p.Results.nfreqs;
specparms.ncyc = p.Results.ncyc;
specparms.space = p.Results.space;
specparms.samplingRate = p.Results.samplingRate;
specparms.showprogress = p.Results.showprogress;
specparms.fvector = p.Results.fvector;
specparms.specnorm = p.Results.specnorm;
specparms.type = p.Results.type;
specparms.winsize = p.Results.winsize;
specparms.overlap = p.Results.overlap;
doPlot = p.Results.doPlot;
saveMat = p.Results.saveMat;
force = p.Results.force;

stringChannels = replace(num2str(lfp.channels),'  ', '-');
stringFreq = replace(num2str(specparms.frange),'  ', '-');
targetFile = dir(['*.power_power_comod_',stringChannels,'_',stringFreq,'Hz.channelinfo.mat']);
if ~isempty(targetFile) && ~force
    disp(['Power power commodulogram for channels ' stringChannels ' in the range ' stringFreq ' Hz already computed. Loading file...']);
    load(targetFile.name);
    return
end

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
    if isfield(lfp,'refCh')
       refCh = lfp.refCh;
    else
        refCh = [];
    end
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

si = 1./samplingRate;


%% Calculate the spectrogram - FFT or WVLT

switch specparms.type
    case 'wavelet'
        %Calcualte the Wavelet Transform
        if size(lfp.data,2) == 1
            [wavespec] = bz_WaveSpec(single(data),...
                'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
            spec = wavespec.data';
        elseif size(lfp.data,2) > 1 
            for i = 1:size(lfp.data,2)
                [wavespec] = bz_WaveSpec(single(data(:,i)),...
                    'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                    'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
                spec{i} = wavespec.data';                
            end
        end
        if ~isempty(refCh)
            [wavespec] = bz_WaveSpec(single(refCh),...
                'frange',specparms.frange,'nfreqs',specparms.nfreqs,'ncyc',specparms.ncyc,...
                'samplingRate',lfp.samplingRate,'space',specparms.space,'fvector',specparms.fvector);
            specRef = wavespec.data';
        else
            specRef = [];
        end
            spectimestamps = timestamps; %Wavelet timestamp are same as LFP        
            comod.freqs = wavespec.freqs;
            
	case 'FFT'
        %Calculate the frequences to use
        if ~isempty(specparms.fvector)
            comod.freqs = fvector;
        else
            switch specparms.space
                case 'log'
                    comod.freqs = logspace(log10(specparms.frange(1)),...
                        log10(specparms.frange(2)),specparms.nfreqs);
                case 'lin'
                    comod.freqs = linspace(specparms.frange(1),...
                        specparms.frange(2),specparms.nfreqs);  
            end
        end
        %Calculate the FFT spectrogram parameters - covert from s to sf
        winsize = specparms.winsize*samplingRate;
        noverlap = specparms.noverlap*samplingRate;
        %Calculate the FFT spectrogram
        if size(lfp.data,2) == 1
            [spec,~,spectimestamps] = spectrogram(single(data),...
                winsize,noverlap,comod.freqs,samplingRate);
        elseif size(lfp.data,2) == 2
            for i = 1:2
                [specT,~,spectimestamps] = spectrogram(single(data(:,i)),...
                    winsize,noverlap,comod.freqs,samplingRate);
                spec{i} = specT;
            end
        end
        if ~isempty(refCh)
            [specRef,~,spectimestamps] = spectrogram(single(refCh),...
                winsize,noverlap,comod.freqs,samplingRate);            
        else
            specRef = [];
        end
            spectimestamps = spectimestamps'+timestamps(1); %account for any time offset
end


%% Calculate the power-power correlations

if isempty(refCh)
    if size(lfp.data,2) == 1
       spec = log10(abs(spec)); %Log-transform power
       comod.corrs = corr(spec','type','spearman');
    elseif size(lfp.data,2) == 2
       for i = 1:2
           spec{i} = log10(abs(spec{i})); %Log-transform power
       end
       comod.corrs = corr(spec{1}',spec{2}','type','spearman');   
    end
elseif exist('refCh')
    specRef = log10(abs(specRef)); %Log-transform power
    if size(lfp.data,2) == 1
       spec = log10(abs(spec)); %Log-transform power
       comod.corrs = corr(specRef',spec','type','spearman');
    elseif size(lfp.data,2) > 1
       for i = 1:size(lfp.data,2)
           spec{i} = log10(abs(spec{i})); %Log-transform power
           comod.corrs{i} = corr(specRef',spec{i}','type','spearman');     
       end
    end    
end


%% Plot
% needs fix for multiple channels

if doPlot && isempty(refCh)  %This whole figure thing can be better.
corrcolor= [makeColorMap([1 1 1],[0 0 0.8],[0 0 0]);...
    makeColorMap([0 0 0],[0.8 0 0],[1 1 1])];
figure
colormap(corrcolor)
imagesc(log2(comod.freqs),log2(comod.freqs),comod.corrs)
% contourf(log2(comod.freqs),log2(comod.freqs),comod.corrs,'LineColor','none');
colorbar
ColorbarWithAxis([-0.4 0.4],'Power-Power Correlation (rho)')
LogScale('xy',2)
xlabel('f (Hz)');ylabel('f (Hz)')
set(gca,'YDir','normal');

powerComod = comod;
powerComod.channels = lfp.channels;
powerComod.specparms = specparms;
    
% NiceSave(['Comodulogram',figparms.plotname],figparms.figfolder,figparms.baseName)
if saveMat
    filename = split(pwd,filesep); filename = filename{end};
    save([filename,'.power_power_comod_',stringChannels,'_',stringFreq,'Hz.channelinfo.mat'],'powerComod');
end

end

