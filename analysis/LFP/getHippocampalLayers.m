
function [hippocampalLayers] = getHippocampalLayers(varargin)
% Identify hippocampal layers based in spectral hallmarks
% 
% INPUTS
% <optional>
% 'basepath'          Default pwd
% 'lfp'               lfp        a buzcode-formatted lfp structure (use bz_GetLFP)
%                       needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%                       If empty or no exist, look for lfp in basePath folder
% 'saveSummary'       Default true
% 'saveMat'           Detault true
% 
% OUTPUT
% hippocampalLayers   channelinfo structure with best channel for stratum
%                       oriens, piramidale, radiatum and slm, and an all
%                       field.
%
% LAYER DEFINITION
% Pyramidal layer is the highest ripple SNR channel.
% Oriens is the highest theta power channel above pyramidal channel.
% Slm is the highest theta power channel below pyramidal channel.
% Radiatum is channel with highest current sink during SPW-ripples.
%
% Manu Valero 2020
% WORK IN PROGRESS!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;

% Deal with inputs
prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.hippocampalLayers.channelinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Hippocampal layers already computed! Loading file.');
    load(targetFile.name);
    return
end

if isempty(lfp)
    try lfp = bz_GetLFP('all');
    catch
        error('No LFP file!!');
    end
end

% Compute channels features
sess = bz_getSessionInfo;
channel_order = sess.AnatGrps.Channels;

% channels.pyramidal = bz_GetBestRippleChan(lfp);
powerProfile_theta = bz_PowerSpectrumProfile([3 12],'channels',[0:31],'showfig',false,'saveMat',false); % [0:63]
powerProfile_hfo = bz_PowerSpectrumProfile([120 250],'channels',[0:31],'showfig',false,'saveMat',false); % [0:63]

channels.slm = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean)));
hfo_theta = powerProfile_hfo.mean - powerProfile_theta.mean;
channels.pyramidal = powerProfile_theta.channels(find(hfo_theta == max(hfo_theta)));

channelsAbovePyr = channel_order(1:find(channel_order==channels.pyramidal)-1);
if length(channelsAbovePyr) == 1
    channels.oriens = NaN;
else
    channels.oriens = 1;% powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean(channelsAbovePyr))));
end

targetFile = dir('*ripples.events.mat');
if ~isempty(targetFile)
    load(targetFile.name);
else
    ripples = bz_FindRipples(bz_GetLFP(channels.pyramidal));
end
twin = 0.1;
[evCsd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
csdRippleProfile = [0 mean(evCsd.data(find(evCsd.timestamps > -10 & evCsd.timestamps < 10),:)) 0];
csdRippleProfile(1:channels.pyramidal) = 0;
channels.radiatum = lfp.channels(find(csdRippleProfile == min(csdRippleProfile)));
channels.all = [channels.oriens; channels.pyramidal; channels.radiatum; channels.slm];

% Summary plot
figure
subplot(1,2,1)
hold on
dev1 = powerProfile_theta.mean - powerProfile_theta.std;
dev2 = powerProfile_theta.mean + powerProfile_theta.std;   
nC = powerProfile_theta.channels;
hold on
fill([nC flip(nC)],[dev1 flip(dev2)],[.8 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
p1 = plot(nC,powerProfile_theta.mean,'color',[.8 .2 .2]);

dev1 = powerProfile_hfo.mean - powerProfile_hfo.std;
dev2 = powerProfile_hfo.mean + powerProfile_hfo.std;
fill([nC flip(nC)],[dev1 flip(dev2)]+ 20,[.2 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
p2 = plot(nC,powerProfile_hfo.mean + 20,'color',[.2 .2 .2]);
xlim([min(nC) max(nC)]);
ax = axis;

plot([channels.pyramidal channels.pyramidal],ax(3:4),'color',[.8 .2 1]);
text(channels.pyramidal, ax(4),'Pyr','color',[.8 .2 1]);

plot([channels.oriens channels.oriens],ax(3:4),'color',[.2 .2 1]);
text(channels.oriens, ax(4),'Or','color',[.2 .2 1]);

plot([channels.radiatum channels.radiatum],ax(3:4),'color',[.5 .5 .1]);
text(channels.radiatum, ax(4),'Rad','color',[.5 .5 .1]);

plot([channels.slm channels.slm],ax(3:4),'color',[.1 .8 .1]);
text(channels.slm, ax(4),'Slm','color',[.1 .8 .1]);

legend([p1 p2], '[3-12 Hz]', '120 240 Hz','Location','southeast');
set(gca,'TickDir','out','XTick',nC); xlabel('Channels'); ylabel('dB');

subplot(1,2,2)
contourf(evCsd.timestamps,(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
box off; colormap(jet); caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);

xs = [evCsd.timestamps(1) evCsd.timestamps(end)];
plot(xs, [channels.pyramidal channels.pyramidal],'color',[.8 .2 1]);
text(xs(2), channels.pyramidal,'Pyr','color',[.8 .2 1]);

plot(xs,[channels.oriens channels.oriens],'color',[.2 .2 1]);
text(xs(2),channels.oriens,'Or','color',[.2 .2 1]);

plot(xs,[channels.radiatum channels.radiatum],'color',[.5 .5 .1]);
text(xs(2),channels.radiatum, 'Rad','color',[.5 .5 .1]);

plot(xs,[channels.slm channels.slm],'color',[.1 .8 .1]);
text(xs(2),channels.slm,'Slm','color',[.1 .8 .1]);
ylim([min(nC) max(nC)]);
set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');

if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,['SummaryFigures\hippocampalLayers.png']);
end

hippocampalLayers = channels;
hippocampalLayers.channels = lfp.channels;

if saveMat
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.hippocampalLayers.channelinfo.mat'],'hippocampalLayers');
end

cd(prevBasepath);
end