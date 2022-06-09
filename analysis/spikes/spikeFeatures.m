function spikeFeatures(varargin)
% Computes and Plot Spike Features (ACG and population CCG)
%
% USAGE
%   spikeFeatures(varargin)
%   
% INPUT (optional)
%   basepath
%   
%   
%
%
% Pablo Abad 2022. Based on computeSessionSummary by Manu Valero 2020
% MV, 2022 exclude and skip intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@bz_isCellInfo);
addParameter(p,'getWaveformsFromDat',false,@islogical);
addParameter(p,'excludeChannels',[],@isnumeric);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
excludeChannels = p.Results.excludeChannels;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;

prevPath = pwd;
cd(basepath);

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat);
end

if skipStimulationPeriods
    try
        optogenetic_responses = getOptogeneticResponse;
        excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
    catch
        warning('Skip stimulation periods not possible...');
    end
end

if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end

%% Spikes Features
disp('Spike-waveform, ACG and cluster location...');

% plot spikes summary
disp('Plotting spikes summary...');
if exist([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'] ,'file')
    load([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords');
    xcoords = chanCoords.x - min(chanCoords.x); xcoords = xcoords/max(xcoords);
    ycoords = chanCoords.y - min(chanCoords.y); ycoords = ycoords/max(ycoords); 
    xcoords(excludeChannels) = [];
    ycoords(excludeChannels) = [];
elseif   exist('chanMap.mat','file')
    load('chanMap.mat','xcoords','ycoords');
    xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
    ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
    xcoords(excludeChannels) = [];
    ycoords(excludeChannels) = [];
else
    xcoords = NaN;
    ycoords = NaN;
end
ccg=CCG(spikes.times,[],'binSize',0.001,'duration',0.08);
dur = 0.08;
xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
figure
set(gcf,'Position',[100 -100 2500 1200])
for jj = 1:size(spikes.UID,2)
    fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
    area(xt,ccg(:,jj,jj),'LineStyle','none');
    try 
        xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,length(spikes.filtWaveform{jj})); % waveform
        spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
        hold on
            plot(xt2,spk,'color','.8 .2 .2');

        plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
        plot((xcoords(spikes.maxWaveformCh1(jj))*30) + dur/2*1000+5+60,...
            ycoords(spikes.maxWaveformCh1(jj))*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
    end
    title(num2str(jj),'FontWeight','normal','FontSize',10);
    if jj == 1
        ylabel('Counts/ norm');
    elseif jj == size(spikes.UID,2)
        xlabel('Time (100 ms /1.5ms)');
    else
        set(gca,'YTick',[],'XTick',[]);
    end
end
saveas(gcf,'SummaryFigures\spikes.png');


cd(prevPath);

end

