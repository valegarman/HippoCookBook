function [phaseMod] = computePhaseModulation(varargin)
%   Calculates distribution of spikes over phases of different oscillations
%   (ripples, SW, theta, low gamma and high gamma)
%
% USAGE
% [rippleModulation, SWModulation, thetaModulation, lgamaModulation, hgammaModulation] = computePhaseModulation
% INPUTS
%   spikes - spike time cellinfo struct
%   rippleChannel - channel for ripple
%   SWChannel - channel for SW, by default radiatum from geHippocampalLayers
%   thetaChannel - channel for theta, by default take oriens from getHippocampallayers
%   hgammaChannel - channel for high gamma, by default oriens from geHippocampalLayers
%   lgammaChannel - channel for low gamma, by default oriens from geHippocampalLayers
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
addParameter(p,'spikes',[]);
addParameter(p,'bandsToCompute',{'rippleModulation','SWModulation','thetaModulation','lgammaModulation',...
    'hgammaModulation','thetaRunModulation','thetaREMModulation'});
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thetaChannel',[],@isnumeric);
addParameter(p,'hgammaChannel',[],@isnumeric);
addParameter(p,'lgammaChannel',[],@isnumeric);
addParameter(p,'ripple_passband',[120 200], @isnumeric);
addParameter(p,'SW_passband',[2 10],@isnumeric);
addParameter(p,'theta_passband',[6 12], @isnumeric);
addParameter(p,'lgamma_passband',[20 60], @isnumeric);
addParameter(p,'hgamma_passband',[60 100],@isnumeric);
addParameter(p,'method','hilbert',@ischar);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'powerThresh',0,@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
bandsToCompute = p.Results.bandsToCompute;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thetaChannel = p.Results.thetaChannel;
hgammaChannel = p.Results.hgammaChannel;
lgammaChannel = p.Results.lgammaChannel;
ripple_passband = p.Results.ripple_passband;
SW_passband = p.Results.SW_passband;
theta_passband = p.Results.theta_passband;
lgamma_passband = p.Results.lgamma_passband;
hgamma_passband = p.Results.hgamma_passband;
method = p.Results.method;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;
powerThresh = p.Results.powerThresh;

%% Session template
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
%% Spikes
if isempty(spikes)
    spikes = loadSpikes;
end

if skipStimulationPeriods
    try
        optogenetic_responses = getOptogeneticResponse;
    catch
        warning('Skip stimulation periods not possible...');
    end
end
excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end

if isempty(powerThresh)
    useThresh = false;
else
    useThresh = true;
end

%% Channels selection
% Ripples
if any([isempty(rippleChannel) isempty(SWChannel) isempty(thetaChannel) isempty(lgammaChannel) isempty(hgammaChannel)])
    [hippocampalLayers] = getHippocampalLayers;
end

if isempty(rippleChannel)
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

if isempty(SWChannel)
    SWChannel = hippocampalLayers.bestShankLayers.radiatum;
end

if isempty(thetaChannel)
    thetaChannel = hippocampalLayers.bestShankLayers.oriens;
end

if isempty(lgammaChannel)
    lgammaChannel = hippocampalLayers.bestShankLayers.oriens;
end

if isempty(hgammaChannel)
    hgammaChannel = hippocampalLayers.bestShankLayers.oriens;
end

%% Initialization of output
rippleMod = [];
SWMod = [];
thetaMod = [];
lgammaMod = [];
hgammaMod = [];
thetaRunMod = [];
thetaREMMod = [];

%% 1. Phase Locking of Spikes to Ripples
if ismember('rippleModulation',bandsToCompute)
    try
        if ~isempty(dir([session.general.name,'.ripples.events.mat']))
            disp('Ripples detected. Loading file !');
            file = dir([session.general.name,'.ripples.events.mat']);
            load(file.name);
        end
        lfpRipple = getLFP(rippleChannel);
        rippleMod = phaseModulation(spikes,lfpRipple,ripple_passband,'intervals',ripples.timestamps,...
            'useThresh',useThresh,'useMinWidth',true,'powerThresh',powerThresh,'method',method);
    catch
        warning('Ripple modulation estimation was not possible...');
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
        SWMod = SWphaseModulation(spikes,lfpSW,SW_passband,'plotting',false,'method',method);
    catch
        warning('Sharp Wave modulation estimation was not possible...');
    end
end

%% 3. Theta Modulation
try thetaEpochs = detectThetaEpochs;
    if ~isfield(thetaEpochs,'thetaREM')
        thetaEpochs = detectThetaEpochs('force',true);
    end
catch
    warning('Detecting theta epochs was not possible. Using full session...');
    thetaEpochs.intervals = [0 Inf];
end


if thetaModulation
    try
        lfpT = getLFP(thetaChannel);
        thetaMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('Sharp Wave modulation estimation was not possible...');
    end
end


%% 4. Low Gamma Modulation
if lgammaModulation
    try
        lfpT = getLFP(lgammaChannel);
        lgammaMod = phaseModulation(spikes,lfpT,lgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
    end
end

%% 5. High Gamma Modulation
if hgammaModulation
    try
        lfpT = getLFP(hgammaChannel);
        hgammaMod = phaseModulation(spikes,lfpT,hgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
    end
end

%% 6. Theta run mod
if thetaRunModulation
    try
        lfpT = getLFP(thetaChannel);
        thetaRunMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.thetaRun.ints,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
    end
end

%% 7. Theta REM mod
if thetaREMModulation
    try
        lfpT = getLFP(thetaChannel);
        thetaREMMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.thetaREM.ints,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
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
    if SWModulation && ~isempty(SWMod)
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

    % ThetaRUN Modulation
    if thetaRunModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([thetaRunMod.phasebins ; thetaRunMod.phasebins + 2*pi],[thetaRunMod.phasedistros(:,i) ;  thetaRunMod.phasedistros(:,i)], 'EdgeColor','none');
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
        saveas(gcf,['SummaryFigures\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end

    % ThetaREM Modulation
    if thetaREMModulation
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([thetaREMMod.phasebins ; thetaREMMod.phasebins + 2*pi],[thetaREMMod.phasedistros(:,i) ;  thetaREMMod.phasedistros(:,i)], 'EdgeColor','none');
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
        saveas(gcf,['SummaryFigures\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end


end


%% Save Output
phaseMod.ripples = rippleMod;
phaseMod.SharpWave = SWMod;
phaseMod.theta = thetaMod;
phaseMod.lgamma = lgammaMod;
phaseMod.hgamma = hgammaMod;
phaseMod.thetaRunMod = thetaRunMod;
phaseMod.thetaREMMod = thetaREMMod;

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
    % ThetaRun
    if thetaRunModulation
        save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'thetaRunMod');
    end
    % thetaREM
    if thetaRunModulation
        save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.PhaseLockingData.cellinfo.mat'],'thetaREMMod');
    end
end

end

