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
% 'force'             Default false
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
% Pablo Abad, Manu Valero 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;

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
    try lfp = getLFP('all');
    catch
        error('No LFP file!!');
    end
end

% Compute channels features
% Updated by Pablo Abad to use session metadata instead of sessionInfo
session = sessionTemplate(basepath,'showGUI',false);
channel_order = session.channels;
% channels.pyramidal = bz_GetBestRippleChan(lfp);
powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'saveMat',false); 
powerProfile_gamma = powerSpectrumProfile(gamma_bandpass,'showfig',true,'saveMat',false);
powerProfile_hfo = powerSpectrumProfile(hfo_bandpass,'showfig',true,'saveMat',false); 

%% Computing Hippocampal Layers by looking at powerSpectrum profiles for theta and hfo
figure
set(gcf,'Position',[100 100 1400 600])
nShanks = length(session.extracellular.spikeGroups.channels);
index = reshape(1:2*nShanks, nShanks, 2).';
for i = 1:length(session.extracellular.spikeGroups.channels)
    if ~all(ismember(session.extracellular.spikeGroups.channels{i},session.channelTags.Bad.channels))
        % Channel slm : Channel with bigger theta amplitude
        channels{i}.slm = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}))))); % 1-index, corresponds to -1 in xml
        hfo_theta = powerProfile_hfo.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i})) - powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        % Channel pyramidal : Channel with bigger amplitude for hfo - theta;
        ch = find(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        channels{i}.pyramidal = powerProfile_theta.channels(ch(find(hfo_theta == max(hfo_theta)))); % 1-index, corresponds to /1 in xml

        channelsAbovePyr = session.extracellular.spikeGroups.channels{i}(1:find(ismember(session.extracellular.spikeGroups.channels{i},channels{i}.pyramidal))-1);
        if isempty(channelsAbovePyr)
            channelsAbovePyr = NaN;
        end
        channels{i}.channelsAbovePyr = channelsAbovePyr;

        if length(channelsAbovePyr) == 1
            channels{i}.oriens = NaN;
        else
            channels{i}.oriens = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean(ismember(powerProfile_theta.channels,sort(channelsAbovePyr))))));
        end
        
        %% Ripples
        ripples{i} = findRipples(channels{i}.pyramidal,'thresholds',[2 5],'passband',[80 200],'durations',[20 150],'saveMat',false);
        
        twin = 0.1;
        [evCsd,lfpAvg] = bz_eventCSD(lfp,ripples{i}.peaks,'channels',session.extracellular.spikeGroups.channels{i},'twin',[twin twin],'plotLFP',false,'plotCSD',false);
        csdRippleProfile = [0 mean(evCsd.data(find(evCsd.timestamps > -10 & evCsd.timestamps < 10),:)) 0];
        shank_channels = session.extracellular.spikeGroups.channels{i};
        csdRippleProfile(find(ismember(shank_channels,shank_channels(1:find(shank_channels == channels{i}.pyramidal))))) = 0;
        channels{i}.radiatum = shank_channels(find(csdRippleProfile == min(csdRippleProfile)));
        % Just in case it detects two radiatum channels
        if length(channels{i}.radiatum) > 1
            channels{i}.radiatum = channels{i}.radiatum(1);
        end
        channels{i}.all = [channels{i}.oriens; channels{i}.pyramidal; channels{i}.radiatum; channels{i}.slm];

        
        %% Summary Plot
        subplot(2,nShanks,index(2*i-1))
        title(['Shank : ', num2str(i)])
        hold on
        % Theta
        [lia,locb] = ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i});
        [B,I] = sort(locb (locb > 0));
        dev1 = powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i})) - powerProfile_theta.ic95(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        dev1 = dev1(I);
        dev2 = powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i})) + powerProfile_theta.ic95(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));   
        dev2 = dev2(I);
%         nC = powerProfile_theta.channels(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        nC = 1:length(dev1);
        hold on
        fill([nC flip(nC)],[dev1 flip(dev2)],[.8 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
        ppmean = powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        ppmean = ppmean(I);
        p1 = plot(nC,ppmean,'color',[.8 .2 .2]);
        % Gamma
        dev1 = powerProfile_gamma.mean(ismember(powerProfile_gamma.channels,session.extracellular.spikeGroups.channels{i})) - powerProfile_gamma.ic95(ismember(powerProfile_gamma.channels,session.extracellular.spikeGroups.channels{i}));
        dev1 = dev1(I);
        dev2 = powerProfile_gamma.mean(ismember(powerProfile_gamma.channels,session.extracellular.spikeGroups.channels{i})) + powerProfile_gamma.ic95(ismember(powerProfile_gamma.channels,session.extracellular.spikeGroups.channels{i}));   
        dev2 = dev2(I);
        fill([nC flip(nC)],[dev1 flip(dev2)] + 10,[.2 .8 .2],'FaceAlpha',.2,'EdgeColor','none');
        ppmean2 = powerProfile_gamma.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        ppmean2 = ppmean2(I);
        p2 = plot(nC, ppmean2 + 10,'color',[.2 .8 .2]);
        xlim([min(nC) max(nC)]);
        ax = axis;
        
        % HFO
        dev1 = powerProfile_hfo.mean(ismember(powerProfile_hfo.channels,session.extracellular.spikeGroups.channels{i})) - powerProfile_hfo.ic95(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        dev1 = dev1(I);
        dev2 = powerProfile_hfo.mean(ismember(powerProfile_hfo.channels,session.extracellular.spikeGroups.channels{i})) + powerProfile_hfo.ic95(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        dev2 = dev2(I);
        fill([nC flip(nC)],[dev1 flip(dev2)] + 20,[.2 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
        ppmean3 = powerProfile_hfo.mean(ismember(powerProfile_hfo.channels,session.extracellular.spikeGroups.channels{i}));
        ppmean3 = ppmean3(I);
        p3 = plot(nC, ppmean3 + 20,'color',[.2 .2 .2]);
        xlim([min(nC) max(nC)]);
        ax = axis;
        set(gca,'XTick',[1:length(ppmean)],'XtickLabels',session.extracellular.spikeGroups.channels{i}); 
        
        if ~isnan(channels{i}.pyramidal)
            plot([find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal) find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal)],ax(3:4),'color',[.8 .2 1]);
            text(find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal), ax(4),'Pyr','color',[.8 .2 1]);
        end
        if ~isnan(channels{i}.oriens)
            plot([find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens) find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens)],ax(3:4),'color',[.2 .2 1]);
            text(find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens), ax(4),'Or','color',[.2 .2 1]);
        end
        if ~isnan(channels{i}.radiatum)
            plot([find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum) find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum)],ax(3:4),'color',[.5 .5 .1]);
            text(find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum), ax(4),'Rad','color',[.5 .5 .1]);
        end
        if ~isnan(channels{i}.slm)
            plot([find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm) find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm)],ax(3:4),'color',[.1 .8 .1]);
            text(find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm), ax(4),'Slm','color',[.1 .8 .1]);
        end

        legend([p1 p2], '[3-12 Hz]', '120 240 Hz','Location','southeast');
        set(gca,'TickDir','out','XTick',nC); xlabel('Channels'); ylabel('dB');
        
        
        % Subplot (1,2,2)
        subplot(2,nShanks,index(2*i))
        title(['Shank : ', num2str(i)])
        contourf(evCsd.timestamps,(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
        box off; colormap(jet); caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);
        xs = [evCsd.timestamps(1) evCsd.timestamps(end)];
%         set(gca,'YTick',[nC(2:end-1)],'YtickLabels',session.extracellular.spikeGroups.channels{i}(nC(2:end-1)));
        
        if ~isnan(channels{i}.pyramidal)
            plot(xs, [find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal) find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal)],'color',[.8 .2 1]);
            text(xs(2), find(session.extracellular.spikeGroups.channels{i} == channels{i}.pyramidal),'Pyr','color',[.8 .2 1]);
        end
        if ~isnan(channels{i}.oriens)
            plot(xs,[find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens) find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens)],'color',[.2 .2 1]);
            text(xs(2),find(session.extracellular.spikeGroups.channels{i} == channels{i}.oriens),'Or','color',[.2 .2 1]);
        end
        if ~isnan(channels{i}.radiatum)
            plot(xs,[find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum) find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum)],'color',[.5 .5 1]);
            text(xs(2),find(session.extracellular.spikeGroups.channels{i} == channels{i}.radiatum),'Rad','color',[.5 .5 1]);
        end
        if ~isnan(channels{i}.slm)
            plot(xs,[find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm) find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm)],'color',[.1 .8 1]);
            text(xs(2),find(session.extracellular.spikeGroups.channels{i} == channels{i}.slm),'Slm','color',[.1 .8 1]);
        end
        ylim([min(nC) max(nC)]);
        set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');     
    end
end
% Save Figure
if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,['SummaryFigures\hippocampalLayers.png']);
end


%% Taking best shank
% First we look how many shanks are in CA1
try CA1Shanks = session.brainRegions.CA1.electrodeGroups;
    numCA1Shanks = length(CA1Shanks);
    disp(['There are: ' , num2str(numCA1Shanks), ' Shanks (#', num2str(CA1Shanks),') in CA1']);
catch
    disp('No brain regions detected! Skiping...');
end

for i = 1:length(session.extracellular.spikeGroups.channels)
    if ~all(ismember(session.extracellular.spikeGroups.channels{i},session.channelTags.Bad.channels))
        % Biggest theta Power
        maxTheta(i) = max(powerProfile_theta.mean((ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}))));
        % Biggest hfo oscillations
        maxHfo(i) = max(powerProfile_hfo.mean((ismember(powerProfile_hfo.channels,session.extracellular.spikeGroups.channels{i}))));
        % Number of Ripples
        numRipples(i) = length(ripples{i}.timestamps);
    end
end

maxTheta_zscore = zscore(maxTheta);
maxHfo_zscore = zscore(maxHfo);
numRipples_zscore = zscore(numRipples);

score = maxTheta_zscore + maxHfo_zscore + numRipples_zscore;
[bestScore, bestShank] = max(score);

% Let's check that best Shank chosen in indeed in CA1
try 
    if ismember(bestShank,session.brainRegions.CA1.electrodeGroups)
        disp('Successfully chosen hippocampal best Shank');
    else
        disp('Best chosen Shank is not in CA1 region..Problem !!');
    end
end

hippocampalLayers.layers = channels;
hippocampalLayers.channels = powerProfile_theta.channels;
hippocampalLayers.bestShank = bestShank;

% definition in best channel
hippocampalLayers.bestShankLayers =  channels{bestShank};

if saveMat
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.hippocampalLayers.channelinfo.mat'],'hippocampalLayers');
end

cd(prevBasepath);
end
