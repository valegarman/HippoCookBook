function spikeFeatures = featureExtraction(varargin)
%
%       [spikeFeatures] = getSpikeFeatures(varargin)
%       
% Extract spike features that then will be used to separate neurons in
% different classes.
% 
% <OPTIONALS>
% basepath
% spikes                spikes struct  
% units ID              type of units to get features from (pyramidal,
%                       narror interneuron, wide interneuron).
% cell_metrics          cell_info struct. Options: 'pyr', 'int', nwint,
%                       wwint
% thetaModulation
% gammaModulation
% 
%
% OUTPUT
% Waveform based measures
%   troughTopeak
%   ab_ratio
% ACG and CCG based metrics
%   acg_tau_rise
% SW phase modulation
% FR per states
%   
% spikeFeatures struct
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[],@bz_isCellInfo);
addParameter(p,'unitsID','all',@isstr);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'lgamma_bandpass',[20 60], @isnumeric);
addParameter(p,'hgamma_bandpass',[60 100], @isnumeric);
addParameter(p,'SW',[],@isstruct);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'SWpassband',[2 10], @isnumeric);
addParameter(p,'ripples',[],@isstruct);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'ripplepassband',[80 200],@isnumeric);
addParameter(p,'plotting',true,@islogical);


parse(p,varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
cell_metrics = p.Results.cell_metrics;
theta_bandpass = p.Results.theta_bandpass;
lgamma_bandpass = p.Results.lgamma_bandpass;
hgamma_bandpass = p.Results.hgamma_bandpass;
unitsID = p.Results.unitsID;
SW = p.Results.SW;
SWChannel = p.Results.SWChannel;
SWpassband = p.Results.SWpassband;
ripples = p.Results.ripples;
rippleChannel = p.Results.rippleChannel;
ripplepassband = p.Results.ripplepassband;
plotting = p.Results.plotting;

keyboard
%% Session Template
session = sessionTemplate(basepath,'showGUI',false);
%% Spikes
disp('Loading Spikes...')
spikes = loadSpikes;

%% Cell_metrics
if  ~isempty(dir([session.general.name,'.cell_metrics.cellinfo.mat']))
    disp('cell_metrics detected ! Loading file.')
    file = dir([session.general.name,'.cell_metrics.cellinfo.mat']);
    load(file.name)
else
    cell_metrics = ProcessCellMetrics('session', session);
end

if ischar(unitsID) && strcmpi(unitsID,'all')
    unitsID = spikes.UID;
elseif ischar(unitsID)
    for i = 1:length(cell_metrics.putativeCellType)
        if strcmpi(cell_metrics.putativeCellType{i},'Pyramidal Cell')
            celltype(i) = 1;
        elseif strcmpi(cell_metrics.putativeCellType{i},'Narrow Interneuron')
            celltype(i) = 2;
        elseif strcmpi(cell_metrics.putativeCellType{i},'Wide Interneuron')
            celltype(i) = 3;
        else
            celltype(i) = 4;
        end
    end
    
    switch unitsID
        case 'pyr'
            unitsID = find(celltype == 1);
        case 'int'
            unitsID = find(celltype == 2 | celltype == 3);
        case 'nwint'
            unitsID = find(celltype == 2);
        case 'wwint'
            unitsID = find(celltype == 3);
    end
end
    
%% Get optogenetic responses to know if a cell is statistically responding to stimulation
try
    if ~isempty(dir([session.general.name,'.optogeneticResponse.cellinfo.mat']))
        disp('Optogenetic Response detected. Loading file !');
        file = dir([session.general.name,'.optogeneticResponse.cellinfo.mat']);
        load(file.name);
    else
        optogeneticResponses = getOptogeneticResponse('numRep',100);
    end
catch
    disp('Not possible to get Optogenetic Responses..');
end
    
%% Phase Locking of Spikes to SharpWaves
% Loading SharpWaves
if isempty(SW)
    if ~isempty(dir([session.general.name,'.sharpwaves.events.mat']))
        disp('sharpwaves detected. Loading file...')
        file = dir([session.general.name,'.sharpwaves.events.mat']);
        load(file.name)
    else
        disp('sharpwaves not found. Running rippleMasterDetector');
        [ripples,SW] = rippleMasterDetector('SWChannel',SWChannel);
    end
end

if isempty(SWChannel) && ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
    file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
    load(file.name);
end
lfpSW = getLFP(SWChannel);
SWMod = SWphaseModulation(spikes,lfpSW,SWpassband,'plotting',false);

if plotting
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    for i = 1:length(spikes.UID)
        subplot(7,ceil(size(spikes.UID,2)/7),i); 
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
    saveas(gcf,['SummaryFigures\SW_',num2str(SWpassband(1)),'-',num2str(SWpassband(end)),'_PhaseModulation.png']);
end

%% Phase Locking of Spikes to Ripples
if isempty(ripples)
    if ~isempty(dir([session.general.name,'.ripples.events.mat']))
        disp('Ripples detected. Loading file !');
        file = dir([session.general.name,'.ripples.events.mat']);
        load(file.name);
    end
end

if isempty(rippleChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers('force',true);
    end
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

lfpRipple = getLFP(rippleChannel);
rippleMod = phaseModulation(spikes,lfpRipple,ripplepassband,'intervals',ripples.timestamps,'useThresh',false,'useMinWidth',false);
if plotting
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    for i = 1:length(spikes.UID)
        subplot(7,ceil(size(spikes.UID,2)/7),i); 
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
    saveas(gcf,['SummaryFigures\ripple_',num2str(ripplepassband(1)),'-',num2str(ripplepassband(end)),'_PhaseModulation.png']);
end

%% Phase Locking to Theta Rhythm
thetaEpochs = detectThetaEpochs;
lfpT = getLFP(thetaEpochs.channel);
thetaMod = phaseModulation(spikes,lfpT,theta_bandpass,'intervals',thetaEpochs.intervals,'useThresh',false,'useMinWidth',false);
if plotting
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    for i = 1:length(spikes.UID)
        subplot(7,ceil(size(spikes.UID,2)/7),i); 
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
    saveas(gcf,['SummaryFigures\theta_',num2str(theta_bandpass(1)),'-',num2str(theta_bandpass(end)),'_PhaseModulation.png']);
end

%% Phase Locking to Low Gamma Rhythm
thetaEpochs = detectThetaEpochs;
lfpT = getLFP(thetaEpochs.channel);
gammaMod = phaseModulation(spikes,lfpT,lgamma_bandpass,'intervals',thetaEpochs.intervals,'useThresh',false,'useMinWidth',false);
if plotting
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    for i = 1:length(spikes.UID)
        subplot(7,ceil(size(spikes.UID,2)/7),i); 
        area([gammaMod.phasebins ; gammaMod.phasebins + 2*pi],[gammaMod.phasedistros(:,i) ;  gammaMod.phasedistros(:,i)], 'EdgeColor','none');
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
    saveas(gcf,['SummaryFigures\lowGamma_',num2str(lgamma_bandpass(1)),'-',num2str(lgamma_bandpass(end)),'_PhaseModulation.png']);
end

%% Phase Locking to High Gamma Rhythm
thetaEpochs = detectThetaEpochs;
lfpT = getLFP(thetaEpochs.channel);
gammaMod = phaseModulation(spikes,lfpT,hgamma_bandpass,'intervals',thetaEpochs.intervals,'useThresh',false,'useMinWidth',false);
if plotting
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    for i = 1:length(spikes.UID)
        subplot(7,ceil(size(spikes.UID,2)/7),i); 
        area([gammaMod.phasebins ; gammaMod.phasebins + 2*pi],[gammaMod.phasedistros(:,i) ;  gammaMod.phasedistros(:,i)], 'EdgeColor','none');
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
    saveas(gcf,['SummaryFigures\highGamma_',num2str(hgamma_bandpass(1)),'-',num2str(hgamma_bandpass(end)),'_PhaseModulation.png']);
end

%% Extraction of Firing Rates in different States
firingRate = [];
firingRate.fr = cell_metrics.firingRate;
firingRate.fr_NREM = cell_metrics.firingRate_NREMstate;
firingRate.fr_REM = cell_metrics.firingRate_REMstate;
firingRate.fr_WAKE = cell_metrics.firingRate_WAKEstate;


%% Waveform-based measures
waveform = [];
waveform.ab_ratio = cell_metrics.ab_rato;
waveform.troughToPeak = cell_metrics.troughToPeak;

%% ACG and CCG based metrics
acg = [];
acg.acg_tau_rise = cell_metrics.acg_tau_rise;
acg.acg_tau_decay = cell_metrics.acg_tau_decay;







%% OUTPUT
spikeFeatures = [];
end

