function [hippocampalLayers] = getHippocampalLayers(varargin)
% [hippocampalLayers] = getHippocampalLayers(varargin)
%
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
% 'promt'             Default false
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

p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'promt',false, @islogical);
addParameter(p,'removeRipplesStimulation',true,@islogical);
addParameter(p,'restrict',[],@isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
theta_bandpass = p.Results.theta_bandpass;
gamma_bandpass = p.Results.gamma_bandpass;
hfo_bandpass = p.Results.hfo_bandpass;
promt = p.Results.promt;
removeRipplesStimulation = p.Results.removeRipplesStimulation;
restrict = p.Results.restrict;

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

if removeRipplesStimulation && ~isempty(dir('*optogeneticPulses.events.mat'))
    targetFile = dir('*optogeneticPulses.events.mat'); load(targetFile.name);
    restrict_temp = SubtractIntervals([0 Inf],optoPulses.stimulationEpochs);
    restrict =  ConsolidateIntervals([restrict; restrict_temp; restrict_temp]);
end

% Compute channels features
% Updated by Pablo Abad to use session metadata instead of sessionInfo
session = loadSession(basepath);
channel_order = session.channels;
% channels.pyramidal = bz_GetBestRippleChan(lfp);
powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'saveMat',true); 
powerProfile_gamma = powerSpectrumProfile(gamma_bandpass,'showfig',true,'saveMat',true);
powerProfile_hfo = powerSpectrumProfile(hfo_bandpass,'showfig',true,'saveMat',true); 

%% Computing Hippocampal Layers by looking at powerSpectrum profiles for theta and hfo
f = figure;
set(gcf,'Position',[100 100 1400 600])
nShanks = length(session.extracellular.spikeGroups.channels);
index = reshape(1:2*nShanks, nShanks, 2).';

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
for i = 1:length(session.extracellular.spikeGroups.channels)
    if ~all(ismember(session.extracellular.spikeGroups.channels{i},session.channelTags.Bad.channels))
        % Channel slm : Channel with bigger theta and gamma amplitude
        theta_gamma = zscor_xnan(powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}))) + zscor_xnan(powerProfile_gamma.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i})));
        ch = find(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        channels{i}.slm = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}))))); % 1-index, corresponds to -1 in xml
        
        % Channel pyramidal : Channel with bigger amplitude for hfo - theta/2;
        hfo_theta = zscor_xnan(powerProfile_hfo.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}))) - zscor_xnan(powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i})));
        ch = find(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
        channels{i}.pyramidal = powerProfile_theta.channels(ch(find(hfo_theta == max(hfo_theta)))); % 1-index, corresponds to /1 in xml
    
        % channels above pyr
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
        
        % if slm is above pyr, then is wrong
        if ismember(channels{i}.slm,channelsAbovePyr)
            channels{i}.slm = NaN;
        end
        
        % Ripples [2 5]
        ripples{i} = findRipples(channels{i}.pyramidal,'thresholds',[2 5],'passband',[80 200],'durations',[20 150],'saveMat',false,'restrict',restrict);
        if isempty(ripples{i}.peaks)
            [~,maxHFOChannel] = max(powerProfile_hfo.mean(session.extracellular.spikeGroups.channels{i}));
            maxHFOChannel = session.extracellular.spikeGroups.channels{i}(maxHFOChannel);
            ripples{i} = findRipples(maxHFOChannel,'thresholds',[.01 .05],'passband',[80 200],'durations',[20 150],'saveMat',false,'restrict',restrict);
        end
        
        twin = 0.1;
        [evCsd,lfpAvg] = bz_eventCSD(lfp,ripples{i}.peaks,'channels',session.extracellular.spikeGroups.channels{i},'twin',[twin twin],'plotLFP',false,'plotCSD',false);
        csdRippleProfile = [NaN mean(evCsd.data(find(evCsd.timestamps > -10 & evCsd.timestamps < 10),:)) NaN];
        shank_channels = session.extracellular.spikeGroups.channels{i};
        csdRippleProfile(find(ismember(shank_channels,shank_channels(1:find(shank_channels == channels{i}.pyramidal))))) = NaN;
        channels{i}.radiatum = shank_channels(find(csdRippleProfile == min(csdRippleProfile)));
        % Just in case it detects two radiatum channels
        if length(channels{i}.radiatum) > 1
            channels{i}.radiatum = channels{i}.radiatum(1);
        end
        if isempty(channels{i}.radiatum)
            channels{i}.radiatum = NaN;
        end

        channels{i}.all = [channels{i}.oriens; channels{i}.pyramidal; channels{i}.radiatum; channels{i}.slm];
        %try
            if promt
                fp = figure;
                set(gcf,'Position',[100 100 500 500]);
                [lia,locb] = ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i});
                [B,I] = sort(locb (locb > 0));
                hold on
                ppmean = powerProfile_theta.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
                ppmean = ppmean(I);
                nC = 1:length(ppmean);
                p1 = plot(ppmean,flip(nC),'color',[.8 .2 .2]);
                plot(ppmean,flip(nC),'o','color',[.8 .2 .2],'MarkerFaceColor',[.8 .2 .2]);
                % Gamma
                ppmean2 = powerProfile_gamma.mean(ismember(powerProfile_theta.channels,session.extracellular.spikeGroups.channels{i}));
                ppmean2 = ppmean2(I);
                p2 = plot(ppmean2 + 10,flip(nC),'color',[.2 .8 .2]);
                plot(ppmean2 + 10,flip(nC),'o','color',[.2 .8 .2],'MarkerFaceColor',[.2 .8 .2]);
                % HFO
                ppmean3 = powerProfile_hfo.mean(ismember(powerProfile_hfo.channels,session.extracellular.spikeGroups.channels{i}));
                ppmean3 = ppmean3(I);
                p3 = plot(ppmean3 + 20,flip(nC),'color',[.2 .2 .2]);
                plot(ppmean3 + 20,flip(nC),'o','color',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2]);
                ylim([min(nC) max(nC) + 4]);
                ax = axis;
                set(gca,'YTick',[1:length(ppmean)],'YtickLabels',flip(session.extracellular.spikeGroups.channels{i})); 

                contourf((evCsd.timestamps  + max(evCsd.timestamps))/14 + ax(2),(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
                targetWin = evCsd.timestamps > -10 & evCsd.timestamps < 10;
                csd_profile = (mean(evCsd.data(targetWin,:)) - nanmean(mean(evCsd.data(targetWin,:))))/nanstd(mean(evCsd.data(targetWin,:)));
                tt = (evCsd.timestamps  + max(evCsd.timestamps))/14 + ax(2); tt = tt(int32(length(tt)/2));
                plot([nan tt + csd_profile nan], nC,'k');
                colormap(jet); caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);
                ax = axis;
                title(['Shank ', num2str(i), ' of ', num2str(length(session.extracellular.spikeGroups.channels))],'FontWeight','normal');

                unsigned_area = [max(nC)+0.5  max(nC)+0.5+4];
                fill([ax([1 2 2 1 1])],[unsigned_area([1 1 2 2 1])],[.7 .7 .7],'EdgeColor','none');
                text(ax(1)+0.5,unsigned_area(2)-1,'unassigned area','FontSize',12,'color',[1 1 1]);


                if ~isnan(channels{i}.pyramidal)
                    l1 = plot(ax(1:2),ones(2,1)*find(flip(session.extracellular.spikeGroups.channels{i}) == channels{i}.pyramidal),'color',[.8 .2 1],'LineWidth',2);
                    t1 = text(l1.XData(2),l1.YData(2),'Pyr','color',[.8 .2 1],'FontSize',12);
                else
                    l1 = plot(ax(1:2),ones(2,1)*max(nC) + 0,'color',[.8 .2 1],'LineWidth',2);
                    t1 = text(l1.XData(2),l1.YData(2),'Pyr','color',[.8 .2 1],'FontSize',12);
                end

                if ~isnan(channels{i}.oriens)
                    l2 = plot(ax(1:2),ones(2,1)*find(flip(session.extracellular.spikeGroups.channels{i}) == channels{i}.oriens),'color',[.2 .2 1],'LineWidth',2);
                    t2 = text(l2.XData(2),l2.YData(2),'Or','color',[.2 .2 1],'FontSize',12);
                else
                    l2 = plot(ax(1:2),ones(2,1)*max(nC) + 1,'color',[.2 .2 1],'LineWidth',2);
                    t2 = text(l2.XData(2),l2.YData(2),'Or','color',[.2 .2 1],'FontSize',12);
                end

                if ~isnan(channels{i}.radiatum)
                    l3 = plot(ax(1:2),ones(2,1)*find(flip(session.extracellular.spikeGroups.channels{i}) == channels{i}.radiatum),'color',[.5 .5 .1],'LineWidth',2);
                    t3 = text(l3.XData(2),l3.YData(2),'Rad','color',[.5 .5 .1],'FontSize',12);
                else
                    l3 = plot(ax(1:2),ones(2,1)* max(nC) + 2,'color',[.5 .5 .1],'LineWidth',2);
                    t3 = text(l3.XData(2),l3.YData(2),'Rad','color',[.5 .5 .1],'FontSize',12);
                end
                if ~isnan(channels{i}.slm)
                    l4 = plot(ax(1:2),ones(2,1)*find(flip(session.extracellular.spikeGroups.channels{i}) == channels{i}.slm),'color',[.8 .1 .1],'LineWidth',2);
                    t4 = text(l4.XData(2),l4.YData(2),'Slm','color',[.8 .1 .1]);
                else
                    l4 = plot(ax(1:2),ones(2,1)* max(nC) + 3,'color',[.8 .1 .1],'LineWidth',2);
                    t4 = text(l4.XData(2),l4.YData(2),'Slm','color',[.8 .1 .1],'FontSize',12);
                end

                draggable(l1,'constraint','v');
                draggable(l2,'constraint','v');
                draggable(l3,'constraint','v');
                draggable(l4,'constraint','v');
                leg = legend([p1 p2 p3], '3-12 Hz', '20-90 Hz','120-240 Hz','Location','southeast');
                set(gca,'TickDir','out','XTick',nC,'XTickLabelRotation',45); ylabel('Channels'); xlabel('dB                                    Time');

                btn = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
                            'Units','normalize','Position', [.82 .93 .15 .06],...
                            'Callback', 'uiresume(gcbf)','BackgroundColor',[.5 .8 .5]);
                disp('Press Accept when done...');
                uiwait(gcf);
                pyr=round(mean(get(l1,'YData')));
                or=round(mean(get(l2,'YData')));
                rad=round(mean(get(l3,'YData')));
                slm=round(mean(get(l4,'YData')));

                delete([l1 l2 l3 l4 t1 t2 t3 t4 leg]);

                l1 = plot(ax(1:2),ones(2,1)*pyr,'color',[.8 .2 1],'LineWidth',2);
                t1 = text(l1.XData(2),l1.YData(2),'Pyr','color',[.8 .2 1],'FontSize',12);

                l2 = plot(ax(1:2),ones(2,1)*or,'color',[.2 .2 1],'LineWidth',2);
                t2 = text(l2.XData(2),l2.YData(2),'Or','color',[.2 .2 1],'FontSize',12);

                l3 = plot(ax(1:2),ones(2,1)*rad,'color',[.5 .5 .1],'LineWidth',2);
                t3 = text(l3.XData(2),l3.YData(2),'Rad','color',[.5 .5 .1],'FontSize',12);

                l4 = plot(ax(1:2),ones(2,1)*slm,'color',[.8 .1 .1],'LineWidth',2);
                t4 = text(l4.XData(2),l4.YData(2),'Slm','color',[.8 .1 .1],'FontSize',12);

                pause(2);
                close(fp);

                channel_list = [flip(session.extracellular.spikeGroups.channels{i}) NaN(1,10)];
                channels{i}.pyramidal = channel_list(pyr)
                channels{i}.radiatum = channel_list(rad);
                channels{i}.slm = channel_list(slm);
                channels{i}.oriens = channel_list(or);
            end
        % end
        
        %% Summary Plot
        subplot(2,nShanks,index(2*i-1))
        title(['Shank : ', num2str(i)]);
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
        set(gca,'TickDir','out','XTick',nC,'XTickLabelRotation',45); xlabel('Channels'); ylabel('dB');
        
        
        % Subplot (1,2,2)
        try
            subplot(2,nShanks,index(2*i))
            title(['Shank : ', num2str(i)]);
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

if promt
    fp = figure;
    hold on
    plot(maxTheta_zscore,'color',[.8 .2 .2],'LineWidth',1.5);
    text(1.1, maxTheta_zscore(1),'Theta score','color',[.8 .2 .2],'FontSize',12);
    plot(maxHfo_zscore,'color',[.3 .3 .3],'LineWidth',1.5);
    text(1.5, maxHfo_zscore(1),'HFO score','color',[.3 .3 .3],'FontSize',12);
    plot(numRipples_zscore,'color',[.5 .5 .9],'LineWidth',1.5);
    text(2, numRipples_zscore(2),'# Ripples score','color',[.5 .5 .9],'FontSize',12);
    
    plot(score,'color',[.4 .9 .4],'LineWidth',2);
    text(1, score(1),'Total score','color',[.4 .9 .4],'FontSize',12);
    
    xlabel('Shanks');     
    
    ax = axis;
    bsl = plot([bestShank bestShank], [ax([3 4])],'color',[.7 1 .7],'LineWidth',1.5);
    bst = text(bestShank+0.1, ax(3)+0.3,'Best shank','color',[.7 1 .7],'FontSize',12);
    
    draggable(bsl,'constraint','h');
    
    btn = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
                        'Units','normalize','Position', [.82 .93 .15 .06],...
                        'Callback', 'uiresume(gcbf)','BackgroundColor',[.5 .8 .5]);
    disp('Press Accept when done...');
    uiwait(gcf);
    
    bestShank=round(mean(get(bsl,'XData')));
    
    delete([bsl bst]);
            
    plot([bestShank bestShank], [ax([3 4])],'color',[.7 1 .7],'LineWidth',1.5);
    text(bestShank+0.1, ax(3)+0.3,'Best shank','color',[.7 1 .7],'FontSize',12);
    
            
    pause(2);
    saveas(gcf,['SummaryFigures\hippocampalLayers_bestShank.png']);
    close(fp);
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
