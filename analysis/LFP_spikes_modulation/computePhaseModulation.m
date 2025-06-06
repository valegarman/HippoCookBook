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
addParameter(p,'bandsToCompute',{'rippleModulation','thetaModulation','lgammaModulation',...
    'hgammaModulation', 'fgammaModulation','thetaRunModulation','thetaREMModulation'});
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thetaChannel',[],@isnumeric);
addParameter(p,'hgammaChannel',[],@isnumeric);
addParameter(p,'lgammaChannel',[],@isnumeric);
addParameter(p,'fgammaChannel',[],@isnumeric);
addParameter(p,'ripple_passband',[120 200], @isnumeric);
addParameter(p,'SW_passband',[2 10],@isnumeric);
addParameter(p,'theta_passband',[6 12], @isnumeric);
addParameter(p,'lgamma_passband',[20 50], @isnumeric);
addParameter(p,'hgamma_passband',[50 100],@isnumeric);
addParameter(p,'fgamma_passband',[100 140],@isnumeric);
addParameter(p,'method','hilbert',@ischar);
addParameter(p,'plotting',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'powerThresh',0,@isnumeric);
addParameter(p,'restrictIntervals',[],@isnumeric); % this is a synonim of restrict_to, kept for funcionality
addParameter(p,'restrict_to',[0 Inf],@isnumeric);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','PhaseLockingData',@ischar);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
bandsToCompute = p.Results.bandsToCompute;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thetaChannel = p.Results.thetaChannel;
hgammaChannel = p.Results.hgammaChannel;
lgammaChannel = p.Results.lgammaChannel;
fgammaChannel = p.Results.fgammaChannel;
ripple_passband = p.Results.ripple_passband;
SW_passband = p.Results.SW_passband;
theta_passband = p.Results.theta_passband;
lgamma_passband = p.Results.lgamma_passband;
hgamma_passband = p.Results.hgamma_passband;
fgamma_passband = p.Results.fgamma_passband;
method = p.Results.method;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
excludeIntervals = p.Results.excludeIntervals;
powerThresh = p.Results.powerThresh;
% restrict_to = p.Results.restrictIntervals;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;

%% Session template
% session = sessionTemplate(basepath,'showGUI',false);
% session = loadSession(basepath);
%% Spikes
if isempty(spikes)
    spikes = loadSpikes;
end

if skipStimulationPeriods
    try
        try 
            targetFile = dir('*optogeneticPulses*');
            optogenetic_responses = importdata(targetFile.name);
        catch
            warning('Could not open optogeneticPulses file... trying to open optogeneticResponses...');
            optogenetic_responses = getOptogeneticResponse;
        end
        excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
    catch
        warning('Skip stimulation periods not possible...');
    end
end

if exist('optogenetic_responses','var') && isfield(optogenetic_responses, 'stimulationEpochs')
    excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
else
    excludeIntervals = [];
end

if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end

session = loadSession;
ints = [];
if isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = '_PhaseLockingData_post';
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif isfield(session,'epochs') && isfield(session.epochs{1},'behavioralParadigm') && restrict_to_baseline
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [0 session.epochs{ii}.startTime];
            warning('Epoch with manipulations found! Restricting analysis to baseline interval!');
        end
    end
    if isempty(ints)
        ints = [0 Inf];
    end
else
    ints = [0 Inf];
end

restrict_ints = IntersectIntervals([ints; restrict_to]);

if any(restrict_ints ~= [0 Inf])
    warning('Restricting analysis for intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},restrict_ints);
        spikes.times{ii} = spikes.times{ii}(status);
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
    try 
        load([basenameFromBasepath(pwd) '.ripples.events.mat']);
        rippleChannel = ripples.detectorinfo.detectionchannel;
    catch
        rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
    end
end

if isempty(SWChannel)
    SWChannel = hippocampalLayers.bestShankLayers.radiatum;
end

if isempty(thetaChannel)
    thetaChannel = hippocampalLayers.bestShankLayers.oriens;
    if isnan(thetaChannel)
        try
            thetaEpochs = detectThetaEpochs();
            thetaChannel = thetaEpochs.channel;
        catch
            disp('Not possible to load good theta channel...');
        end
    end
end

if isempty(lgammaChannel)
    lgammaChannel = hippocampalLayers.bestShankLayers.radiatum;
    if isnan(lgammaChannel)
        lgammaChannel = rippleChannel;
    end
end

if isempty(hgammaChannel)
    hgammaChannel = hippocampalLayers.bestShankLayers.oriens;
    if isnan(hgammaChannel)
        thetaEpochs = detectThetaEpochs();
        hgammaChannel = thetaEpochs.channel;
    end
end

if isempty(fgammaChannel)
    fgammaChannel = rippleChannel;
    if isnan(fgammaChannel)
        fgammaChannel = hippocampalLayers.bestShankLayers.pyramidal;
    end
end

%% Initialization of output
rippleMod = [];
SWMod = [];
thetaMod = [];
lgammaMod = [];
hgammaMod = [];
fgammaMod = [];
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
if ismember('SWModulation',bandsToCompute)
    try
        if ~isempty(dir([session.general.name,'.sharpwaves.events.mat']))
            disp('Sharpwaves detected. Loading file...');
            file = dir([session.general.name,'.sharpwaves.events.mat']);
            load(file.name)
        end
        
        lfpSW = getLFP(SWChannel);
        SWMod = SWphaseModulation(spikes,lfpSW,SW_passband,'plotting',false);
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

if ismember('thetaModulation',bandsToCompute)
    try
        lfpT = getLFP(thetaChannel);
        thetaMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('Theta modulation estimation was not possible...');
    end
end


%% 4. Low Gamma Modulation
if ismember('lgammaModulation',bandsToCompute)
    try
        lfpT = getLFP(lgammaChannel);
        lgammaMod = phaseModulation(spikes,lfpT,lgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('Low gamma modulation estimation was not possible...');
    end
end

%% 5. High Gamma Modulation
if ismember('hgammaModulation',bandsToCompute)
    try
        lfpT = getLFP(hgammaChannel);
        hgammaMod = phaseModulation(spikes,lfpT,hgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('High gamma modulation estimation was not possible...');
    end
end

%% 6. Fast Gamma Modulation
if ismember('fgammaModulation',bandsToCompute)
    try
        lfpT = getLFP(fgammaChannel);
        fgammaMod = phaseModulation(spikes,lfpT,fgamma_passband,'intervals',thetaEpochs.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('High gamma modulation estimation was not possible...');
    end
end


%% 6. Theta run mod
if ismember('thetaRunModulation',bandsToCompute)
    try
        lfpT = getLFP(thetaChannel);
        thetaRunMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.thetaRun.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('Theta run modulation estimation was not possible...');
    end
end

%% 7. Theta REM mod
if ismember('thetaREMModulation',bandsToCompute)
    try
        lfpT = getLFP(thetaChannel);
        thetaREMMod = phaseModulation(spikes,lfpT,theta_passband,'intervals',thetaEpochs.thetaREM.intervals,...
            'useThresh',useThresh,'useMinWidth',false,'powerThresh',powerThresh,'method',method);
    catch
        warning('Theta REM modulation estimation was not possible...');
    end
end

%% Plotting
if plotting
    % Ripple Modulation
    try
    if ismember('rippleModulation',bandsToCompute)
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
        saveas(gcf,['SummaryFigures\ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),save_as, '.png']);
    end
    
    % Sharpwave Modulation
    if ismember('SWModulation',bandsToCompute) && ~isempty(SWMod)
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
        saveas(gcf,['SummaryFigures\SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)), save_as, '.png']);
    end
    
    % Theta Modulation
    if ismember('thetaModulation',bandsToCompute)
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)), save_as, '.png']);
    end
    
    % LowGama Modulation
    if ismember('lgammaModulation',bandsToCompute)
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(lgammaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\lGamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)), save_as, '.png']);
    end
    
    % High Modulation
    if ismember('hgammaModulation',bandsToCompute)
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(hgammaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\hGamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),save_as ,'.png']);
    end

    % Fast Modulation
    if ismember('fgammaModulation',bandsToCompute)
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for i = 1:length(spikes.UID)
            subplot(7,ceil(size(spikes.UID,2)/7),i)
            area([fgammaMod.phasebins ; fgammaMod.phasebins + 2*pi],[fgammaMod.phasedistros(:,i) ;  fgammaMod.phasedistros(:,i)], 'EdgeColor','none');
            hold on;
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(i),'FontWeight','normal','FontSize',10);
            if i == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(fgammaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\fGamma_',num2str(fgamma_passband(1)),'-',num2str(fgamma_passband(end)),save_as ,'.png']);
    end

    % ThetaRUN Modulation
    if ismember('thetaRunModulation',bandsToCompute)
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)), save_as, '.png']);
    end

    % ThetaREM Modulation
    if ismember('thetaREMModulation',bandsToCompute) && ~isempty(thetaREMMod)
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif i == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,['SummaryFigures\thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)), save_as,'.png']);
    end
    end

end


%% Save Output
phaseMod.ripples = rippleMod;
phaseMod.SharpWave = SWMod;
phaseMod.theta = thetaMod;
phaseMod.lgamma = lgammaMod;
phaseMod.hgamma = hgammaMod;
phaseMod.fgamma = fgammaMod;
phaseMod.thetaRunMod = thetaRunMod;
phaseMod.thetaREMMod = thetaREMMod;

if saveMat
    % Ripple 
    if ismember('rippleModulation',bandsToCompute)
        save([session.general.name,'.ripple_',num2str(ripple_passband(1)),'-',num2str(ripple_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'rippleMod');
    end
    % Sharp Wave
    if ismember('SWModulation',bandsToCompute)
        save([session.general.name,'.SW_',num2str(SW_passband(1)),'-',num2str(SW_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'SWMod');
    end
    % Theta
    if ismember('thetaModulation',bandsToCompute)
        save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'thetaMod');
    end
    % Low Gamma
    if ismember('lgammaModulation',bandsToCompute)
        save([session.general.name,'.lgamma_',num2str(lgamma_passband(1)),'-',num2str(lgamma_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'lgammaMod');
    end
    % High Gamma
    if ismember('hgammaModulation',bandsToCompute)
        save([session.general.name,'.hgamma_',num2str(hgamma_passband(1)),'-',num2str(hgamma_passband(end)),...
            '.', save_as,'.cellinfo.mat'],'hgammaMod');
    end
    % Fast Gamma
    if ismember('fgammaModulation',bandsToCompute)
        save([session.general.name,'.fgamma_',num2str(fgamma_passband(1)),'-',num2str(fgamma_passband(end)),...
            '.', save_as,'.cellinfo.mat'],'fgammaMod');
    end
    % ThetaRun
    if ismember('thetaRunModulation',bandsToCompute)
        save([session.general.name,'.thetaRun_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'thetaRunMod');
    end
    % thetaREM
    if ismember('thetaRunModulation',bandsToCompute)
        save([session.general.name,'.thetaREM_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
            '.', save_as, '.cellinfo.mat'],'thetaREMMod');
    end
end

end

