function [powerProfile] = powerSpectrumProfile(frange,varargin)
% power distribution per channels
%
%   INPUTS
%   frange      (Hz, If empty take by default 6 12, theta)
%
%   (optional)
%    lfp        a buzcode-formatted lfp structure (getLFP)
%               needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%               If empty or no exist, look for lfp in basePath folder
%    winsize    size of the silding time window (s, default 2)
%    dt         sliding time interval (s, default 1)
%    channels   subset of channels to calculate PowerSpectrumSlope
%               (default: all)
%    showfig    true/false - show a summary figure of the results
%               (default:false)
%    saveMat    put your basePath here to save/load
%               baseName.PowerSpectrumProfile_'frange'.lfp.mat  (default: true)
%   forceDetect (default false)
%
%   OUTPUTS
%   powerProfile
%       .mean       log10-transformed mean amplitude of the spectrogram
%       .std
%       .median
%       .channels
%       .channels_shank
%       .frange
%       .f
%       .ic95
%
% MV-BuzsakiLab 2019
% Edited by Peter Petersen
% dited by Pablo Abad (compute NaN on bad channels from sessionTemplate and
%                           removed sessionInfo dependencies)
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Parsing inputs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'winSize',4,@isscalar)
addParameter(p,'dt',2, @isscalar)
addParameter(p,'channels','all')
addParameter(p,'showfig',true,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'lfp',[])
addParameter(p,'forceDetect',false,@islogical)
addParameter(p,'useParfor',true,@islogical)
addParameter(p,'rejectChannels',[],@isnumeric);

parse(p,varargin{:})
showfig = p.Results.showfig;
saveMat = p.Results.saveMat;
channels = p.Results.channels;
winSize = p.Results.winSize;
lfp = p.Results.lfp;
dt = p.Results.dt;
forceDetect = p.Results.forceDetect;
useParfor = p.Results.useParfor;
rejectChannels = p.Results.rejectChannels;
basepath = p.Results.basepath;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Resolving inputs   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
try [session] = loadSession(basepath);
    if (exist([session.general.name,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat'],'file') || ...
            exist([session.general.name,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.lfp.mat'],'file')) ...
            && ~forceDetect
        disp(['Power spectrum profile already calculated for ', session.general.name, ' in the range: ', num2str(frange(1)),' - ',num2str(frange(2)), 'Hz. Loading file.']);
        try load([session.general.name,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat']);
        catch
            load([session.general.name,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.lfp.mat']);
        end
        return
    end
catch disp('No sessionTemplate file!!')
    if ischar('channels') && strcmpi(channels,'all')
        lfp = getLFP('all');
        channels = lfp.channels+1;
    end
    keyboard;
    session.extracellular.srLfp = lfp.samplingRate;
    useParfor = false;
    fileName = strsplit(basepath,filesep);
    fileName = fileName{end};
    session.general.name = fileName;
    session.extracellular.spikeGroups.channels = lfp.channels+1; 
end

if ~exist('frange') || isempty(frange)
    frange = [6 12];
end
  
%% Dealing with channels input
channels = 1:session.extracellular.nChannels;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Calculate spectrogram per channel
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
params.Fs = session.extracellular.srLfp; params.fpass = frange; 
params.tapers = [3 5]; params.pad = 1;
tic
powerProfileMean = [];
powerProfileStd = [];
powerProfileIc95 = [];
powerProfileMedian = [];
powerProfileChannels = [];
disp('Calculating spectrograms channelwise')
if useParfor
    parfor (ii = 1:length(channels),18)
        fprintf('Channel %3.i/%3.i, ',ii, length(channels));
        lfp = getLFP(channels(ii),'noPrompts', true); % now channels are 1-index an getLFP_temp is also 1-index
        [S,t,f] = mtspecgramc_fast(single(lfp.data),[4 2],params);
        S = 10 * log10(S);
        S(find(isinf(S))) = NaN;
        t(find(isinf(t))) = NaN;
        powerProfileMean(ii) = nanmean(nanmean(S,2));
        powerProfileStd(ii) = nanstd(nanmean(S,2));
        powerProfileIc95(ii) = 1.96 * nanstd(nanmean(S,2))/sqrt(length(nanmean(S,2)));
        powerProfileMedian(ii) = nanmedian(nanmedian(S,2));
        powerProfileChannels(ii) = channels(ii);
    end
else
    if isempty(lfp)
        lfp = getLFP('all','noPrompts', true);
    end
    for ii = 1:length(channels)
        fprintf('Channel %3.i/%3.i ,\n',ii, length(channels));
        [S,t,f] = mtspecgramc_fast(single(lfp.data(:,channels(ii))),[4 2],params);
        S = 10 * log10(S);
        S(find(isinf(S))) = NaN;
        t(find(isinf(S))) = NaN;
        powerProfileMean(ii) = nanmean(nanmean(S,2));
        powerProfileStd(ii) = std(nanmean(S,2));
        powerProfileIc95(ii) = 1.96 * nanstd(nanmean(S,2))/sqrt(length(nanmean(S,2)));
        powerProfileMedian(ii) = nanmedian(nanmedian(S,2));
        powerProfileChannels(ii) = channels(ii);
    end
end
toc

% In this part we compute NaNs in bad channels.
% Also, in all the channels that belong to a shank are considered bad,
% remove those channels
powerProfileMean(session.channelTags.Bad.channels) = NaN;
powerProfileStd(session.channelTags.Bad.channels) = NaN;
powerProfileIc95(session.channelTags.Bad.channels) = NaN;
powerProfileMedian(session.channelTags.Bad.channels) = NaN;

powerProfile.mean = powerProfileMean;
powerProfile.std = powerProfileStd;
powerProfile.ic95 = powerProfileIc95;
powerProfile.median = powerProfileMedian;
powerProfile.channels = powerProfileChannels;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Saving the result to basename.PowerSpectrumProfile_frange.channelinfo.mat
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
powerProfile.processinginfo.function = 'bz_PowerSpectrumProfile';
powerProfile.processinginfo.date = now;
powerProfile.processinginfo.params.winSize = winSize;
powerProfile.processinginfo.params.dt = dt;
powerProfile.processinginfo.params.frange = frange;
if saveMat
    save([session.general.name,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat'],'powerProfile');
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if showfig
    figure('Position', get(0,'screensize'),'Name',session.general.name)
    cmap = jet(size(session.extracellular.electrodeGroups.channels,2));
    plt1 = [];
    for ii = 1:size(session.extracellular.electrodeGroups.channels,2)
        if ~all(ismember(session.extracellular.electrodeGroups.channels{ii},session.channelTags.Bad.channels))
            [Lia] = ismember(session.extracellular.electrodeGroups.channels{ii}, channels);
            nC = 1:length(session.extracellular.electrodeGroups.channels{ii});
            nC = nC(Lia);
            hold on

            dev1 = powerProfile.mean(session.extracellular.electrodeGroups.channels{ii}(Lia)) - powerProfile.ic95(session.extracellular.electrodeGroups.channels{ii}(Lia));
            dev2 = powerProfile.mean(session.extracellular.electrodeGroups.channels{ii}(Lia)) + powerProfile.ic95(session.extracellular.electrodeGroups.channels{ii}(Lia));

            hold on
            fill([dev1 flip(dev2)],[nC flip(nC)],cmap(ii,:),'FaceAlpha',.2,'EdgeColor','none')
            try plt1(ii) = plot(powerProfile.mean(session.extracellular.electrodeGroups.channels{ii}(Lia)),nC(Lia),'color',cmap(ii,:)); end
        end
    end
    
    ax=axis; axis tight; xlim(ax(1:2));
    ylabel('Channels'); xlabel('Power'); title(strcat('Freq range:',num2str(frange),'Hz'),'FontWeight','normal');
    set(gca,'YDir','reverse');
%     legend(plt1,{num2str([1:size(session.extracellular.electrodeGroups.channels,2)]')})
    if ~exist('SummaryFigures','dir')
        mkdir('SummaryFigures')
    end
    saveas(gcf,['SummaryFigures\PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.png']);

end

