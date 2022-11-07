function [phaseMod] = computePhaseModulationPerSubSession(varargin)
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
addParameter(p,'includeIntervals',[],@isnumeric);

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
includeIntervals = p.Results.includeIntervals;
%% Session template
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);
%% Spikes
if isempty(spikes)
    spikes = loadSpikes;
end

%% Get MergePoints
try
    targetFile = dir('*MergePoints.events.mat');
    load(targetFile.name);
catch
    error('MergePoints not available. Quitting analysis...');
end

for ii = 1:length(MergePoints.foldernames)
    ts = MergePoints.timestamps(ii,:);
    phaseModulation = computePhaseModulation('rippleChannel',rippleChannel,'SWChannel',SWChannel,'restrictIntervals',ts,'saveMat',false,'plotting',false);
    
    phaseMod.(MergePoints.foldernames{ii}) = phaseModulation;
end


%% Output

if saveMat
    save([session.general.name,'.PhaseLockingDataSubSessions.cellinfo.mat'],'phaseMod','-v7.3');    
end



%% Plotting
flds = fields(phaseMod);

for ii = 1:length(flds)
    % Ripple modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).ripples.phasebins ; phaseMod.(flds{ii}).ripples.phasebins + 2*pi],[phaseMod.(flds{ii}).ripples.phasedistros(:,jj) ;  phaseMod.(flds{ii}).ripples.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(rippleChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        mkdir('SubSessions');
        saveas(gcf,['SubSessions\ripple_',flds{ii},'_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),'_PhaseModulation.png']);
    end
    
    % Sharpwave Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).SharpWave.phasebins ; phaseMod.(flds{ii}).SharpWave.phasebins + 2*pi],[phaseMod.(flds{ii}).SharpWave.phasedistros(:,jj) ;  phaseMod.(flds{ii}).SharpWave.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(SWChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\SW_',flds{ii},'_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),'_PhaseModulation.png']);
    end
    % Theta Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).theta.phasebins ; phaseMod.(flds{ii}).theta.phasebins + 2*pi],[phaseMod.(flds{ii}).theta.phasedistros(:,jj) ;  phaseMod.(flds{ii}).theta.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(phaseMod.(flds{ii}).theta.detectorParams.channels)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\theta_',flds{ii},'_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end
    
    % LowGamma Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).lgamma.phasebins ; phaseMod.(flds{ii}).lgamma.phasebins + 2*pi],[phaseMod.(flds{ii}).lgamma.phasedistros(:,jj) ;  phaseMod.(flds{ii}).lgamma.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(phaseMod.(flds{ii}).theta.detectorParams.channels)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\lGamma_',flds{ii},'_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),'_PhaseModulation.png']);
    end
    
    % HighGamma Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).hgamma.phasebins ; phaseMod.(flds{ii}).hgamma.phasebins + 2*pi],[phaseMod.(flds{ii}).hgamma.phasedistros(:,jj) ;  phaseMod.(flds{ii}).hgamma.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(phaseMod.(flds{ii}).theta.detectorParams.channels)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\hGamma_',flds{ii},'_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),'_PhaseModulation.png']);
    end
    % ThetaRUN Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).thetaRunMod.phasebins ; phaseMod.(flds{ii}).thetaRunMod.phasebins + 2*pi],[phaseMod.(flds{ii}).thetaRunMod.phasedistros(:,jj) ;  phaseMod.(flds{ii}).thetaRunMod.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(phaseMod.(flds{ii}).theta.detectorParams.channels)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\thetaRun_',flds{ii},'_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end
    
    % ThetaREM Modulation
    try
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),jj)
            area([phaseMod.(flds{ii}).thetaREMMod.phasebins ; phaseMod.(flds{ii}).thetaREMMod.phasebins + 2*pi],[phaseMod.(flds{ii}).thetaREMMod.phasedistros(:,jj) ;  phaseMod.(flds{ii}).thetaREMMod.phasedistros(:,jj)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(phaseMod.(flds{ii}).theta.detectorParams.channels)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SubSessions\thetaREM_',flds{ii},'_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),'_PhaseModulation.png']);
    end
end
close all;






end
