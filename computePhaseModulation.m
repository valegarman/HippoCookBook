function [rippleMod, SWMod, thetaMod, lgammaMod, hgammaMod] = computePhaseModulation(varargin)
%   Calculates distribution of spikes over phases of different oscillations
%   (ripples, SW, theta, low gamma and high gamma)
%
% USAGE
% [rippleModulation, SWModulation, thetaModulation, lgamaModulation, hgammaModulation] = computePhaseModulation
% INPUTS
%   spikes - spike time cellinfo struct
%   rippleModulation - logical if want to compute ripple modulation (default true)
%   SWModulation - logical if want compute sharp wave modulation (default true)
%   thetaModulation - logical if want compute theta modulation (default true)
%   lgammaModulation - logical if want compute low gamma Modulation (default true)
%   hgammaModulation - logical if want compute high gamma Modulation (default true)
%   rippleChannel - channel for ripple
%   SWChannel - channel for SW
%   thetaChannel - channel for theta
%   ripple_passband - frequency range for ripple phase modulation
%   SW_passband - frequency range for sharp wave modulation
%   theta_passband - frequency range for theta modulation
%   lgamma_passband - frequency range for low gamma modulation
%   hgamma_passband - frequency range for high gamma modulation
%   plotting - logical if you want to plot, false if not (default true)
%   method - method selection for how to generate phase: 'hilbert'(default),
%       or 'wavelet'
%   saveMat - logical to save cellinfo.mat file with results for each
%       modulation (default true)
%
% OUTPUTS
% phasedistros  - Spike distribution perecentages for each cell in each bin
%               specified by phasebins
%
% phasebins     - 180 bins spanning from 0 to 2pi
%
% phasestats    - ncellsx1 structure array with following (via
%                 CircularDistribution.m from FMAToolbox)
%                    phasestats.m        mean angle
%                    phasestats.mode     distribution mode
%                    phasestats.k        concentration
%                    phasestats.p        p-value for Rayleigh test
%                    phasestats.r        mean resultant length
% Pablo Abad and Manuel Valero 2022
%
%% Defaults and params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@bz_isCellInfo);
addParameter(p,'lfp',[],@bz_isLFP);
addParameter(p,'rippleModulation',true,@islogical);
addParameter(p,'SWModulation',true,@islogical);
addParameter(p,'thetaModulation',true,@islogical);
addParameter(p,'lgammaModulation',true,@islogical);
addParameter(p,'hgammaModulation',true,@islogical);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thetaChannel',[],@isnumeric);
addParameter(p,'ripple_passband',[120 200], @isnumeric);
addParameter(p,'SW_passband',[2 10],@isnumeric);
addParameter(p,'theta_passband',[6 12], @isnumeric);
addParameter(p,'lgamma_passband',[20 60], @isnumeric);
addParameter(p,'hgamma_passband',[60 100],@isnumeric);
addParameter(p,'method','wavelet',@isstr);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
lfp = p.Results.lfp;
rippleModulation = p.Results.rippleModulation;
SWModulation = p.Results.SWModulation;
thetaModulation = p.Results.thetaModulation;
lgammaModulation = p.Results.lgammaModulation;
hgammaModulation = p.Results.hgammaModulation;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thetaChannel = p.Results.thetaChannel;
ripple_passband = p.Results.ripple_passband;
SW_passband = p.Results.SW_passband;
theta_passband = p.Results.theta_passband;
lgamma_passband = p.Results.lgamma_passband;
hgamma_passband = p.Results.hgamma_passband;
method = p.Results.method;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;

%% Session template
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
%% Spikes
if isempty(spikes)
    spikes = loadSpikes;
end

%% Channels selection
% Ripples
if isempty(rippleChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        disp('hippocampal layers found. Loading file...');
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers('force',true);
    end
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

% Sharp Waves
if isempty(SWChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        disp('hippocampal layers found. Loading file...');
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers('force',true);
    end
    SWChannel = hippocampalLayers.bestShankLayers.radiatum;
end
%% Initialization of output
rippleMod = [];
SWMod = [];
thetaMod = [];
lgammaMod = [];
hgammaMod = [];

%% 1. Phase Locking of Spikes to Ripples
if rippleModulation
    try
        if ~isempty(dir([session.general.name,'.ripples.events.mat']))
            disp('Ripples detected. Loading file !');
            file = dir([session.general.name,'.ripples.events.mat']);
            load(file.name);
        end
        
        lfpRipple = getLFP(rippleChannel);
        rippleMod = phaseModulation(spikes,lfpRipple,ripple_passband,'intervals',ripples.timestamps,...
            'useThresh',false,'useMinWidth',false);
    catch
    end  
end


%% 2. Sharp Wave Modulation
if SWModulation
    try
        if ~isempty(dir([session.general.name,'.sharpwaves.events.mat']))
            disp('Sharpwaves detected. Loading file...');
            file = dir([session.general.name,'.sharpwaves.events.mat']);
            load(file.name)
        end
        
        lfpSW = getLFP(SWChannel);
        SWMod = SWphaseModulation(spikes,lfpSW,SW_passband,'plotting',false);
    catch
    end
end

%% 3. Theta Modulation
if thetaModulation
    try
        thetaEpochs = detectThetaEpochs;
        lfpT = getLFP(thetaEpochs.channel);
        thetaMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',false,'useMinWidth',false);
    catch
    end
end


%% 4. Low Gamma Modulation
if lgammaModulation
    try
        thetaEpochs = detectThetaEpochs;
        lfpT = getLFP(thetaEpochs.channel);
        lgammaMod = phaseModulation(spikes,lfpT,lgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',false,'useMinWidth',false);
    catch
    end
end

%% 5. High Gamma Modulation
if hgammaModulation
    try
        thetaEpochs = detectThetaEpochs;
        lfpT = getLFP(thetaEpochs.channel);
        hgammaMod = phaseModulation(spikes,lfpT,hgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',false,'useMinWidth',false);
    catch
    end
end

%% Plotting
if plotting
    % Ripple Modulation
    if rippleModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([rippleMod.phasebins ; rippleMod.phasebins + 2*pi],[rippleMod.phasedistros(:,i) ;  rippleMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_PhaseModulation.png']);
    end
    
    % Sharpwave Modulation
    if SWModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([SWMod.phasebins ; SWMod.phasebins + 2*pi],[SWMod.phasedistros(:,i) ;  SWMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_PhaseModulation.png']);
    end
    
    % Theta Modulation
    if thetaModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([thetaMod.phasebins ; thetaMod.phasebins + 2*pi],[thetaMod.phasedistros(:,i) ;  thetaMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end
    
    % LowGama Modulation
    if lgammaModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([lgammaMod.phasebins ; lgammaMod.phasebins + 2*pi],[lgammaMod.phasedistros(:,i) ;  lgammaMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PhaseModulation.png']);
    end
    
    % High Modulation
    if hgammaModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([hgammaMod.phasebins ; hgammaMod.phasebins + 2*pi],[hgammaMod.phasedistros(:,i) ;  hgammaMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaEpochs.channel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PhaseModulation.png']);
    end
end


%% Save Output

if saveMat
    % Ripple 
    if rippleModulation
        save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'rippleMod');
    end
    % Sharp Wave
    if SWModulation
        save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'SWMod');
    end
    % Theta
    if thetaModulation
        save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'thetaMod');
    end
    % Low Gamma
    if lgammaModulation
        save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'lgammaMod');
    end
    % High Gamma
    if hgammaModulation
        save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'hgammaMod');
    end
end

end

