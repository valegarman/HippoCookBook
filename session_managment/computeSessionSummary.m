
function computeSessionSummary(varargin)
% Runs a list of preliminary descriptive analysis 
%
% USAGE 
%   computeSessionSummary(varargin)
%
% INPUT (optional)
%   basepath         By default pwd
%   listOfAnalysis   Summary analysis that will be computed (as a cell with
%                        strings). Default 'all'. Possible values: 'spikes',
%                        'analogPulses', 'digitalPulses',
%                        'downStates','ripples', 'tMazeBehaviour',
%                        'linearMazeBehaviour', 'thetaModulation'
%   exclude         Cell with strings with the list of analysis to exclude
%  
%   (specific analysis optionss)
%   excludeShanks           Default []
%   analogChannelsList      Array of channel to perform 'analogPulses' psth and csd. Default 'all'
%   digitalChannelsList     Array of channel to perform 'digitalPulses' psth and csd. Default 'all'
%   skipErrors      Default, true
%
% Manu Valero-BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'listOfAnalysis','all');
addParameter(p,'exclude',[]);
addParameter(p,'excludeShanks',[],@isnumeric);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;
tracking_pixel_cm = p.Results.tracking_pixel_cm;

prevPath = pwd;
cd(basepath);

if ischar(listOfAnalysis) && strcmpi(listOfAnalysis,'all')
    listOfAnalysis = {'spikes', 'analogPulses', 'digitalPulses', 'hippocampalLayers','downStates', 'ripples', 'tMazeBehaviour','linearMazeBehaviour','thetaModulation'};
end
if ~isempty(exclude)
    listOfAnalysis(ismember(listOfAnalysis, exclude)) = [];
end
session = loadSession;
session = sessionTemplate(session);
session.channels = 1:session.extracellular.nChannels;
save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');

excludeChannels = [];
for ii = 1:length(excludeShanks)
    excludeChannels = session.extracellular.electrodeGroups.channels{excludeShanks(ii)};
end

mkdir('SummaryFigures'); % create folder
close all

try
    getuLEDsPulses();
end
cd(basepath);
% SPIKES SUMMARY
if any(ismember(listOfAnalysis,'spikes'))
    try
       spikes = loadSpikes('forceReload',false);
       spikeFeatures();
       getAverageCCG;
       clear spikes;
    catch
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end

% ANALOG CHANNELS CSD 
if any(ismember(listOfAnalysis,'analogPulses'))
    try 
        disp('Psth and CSD from analog-in inputs...');
        pulses = getAnalogPulses('manualThr',true);
        if ~isfield(pulses,'analogChannelsList') && ~isempty(pulses)
            pulses.analogChannelsList = ones(size(pulses.amplitude))*65;
            pulses.eventGroupID = ones(size(pulses.amplitude));
            save([basenameFromBasepath '.pulses.events.mat'],'pulses')
        end

        if ~isempty(pulses)
            if ischar(analogChannelsList) && strcmpi(analogChannelsList,'all')
                listOfChannel = unique(pulses.analogChannel);
            else
                listOfChannel = analogChannelsList;
            end
            
            % CSD
            for mm = 1:length(listOfChannel)
                fprintf('Stimulus %3.i of %3.i \n',mm, length(listOfChannel)); %\n
                st = pulses.timestamps(pulses.analogChannelsList == listOfChannel(mm),1);
                if isempty(st)
                    st = pulses.timestamps(pulses.analogChannelsList == listOfChannel(mm) + 1,1);
                end
                shanks = session.extracellular.electrodeGroups.channels;            
                shanks(excludeShanks) = [];
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:length(shanks)
                    lfp = getLFP(shanks{jj});
                    twin = 0.02;
                    [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('AnalogCh #',num2str(listOfChannel(mm))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                saveas(gcf,['SummaryFigures\analogPulsesCSD_ch',num2str(listOfChannel(mm)), '.png']);
                clear pulses listOfChannel
            end
        else
            warning('Analog pulses not found!');
        end

    catch
        warning('Error on CSD from analog inputs! ');
    end
end

% DIGITAL CHANNELS CSD
if any(ismember(listOfAnalysis,'digitalPulses'))
    try 
        disp('Psth and CSD from digital-in inputs...');
        pulses = getDigitalIn;
        if ~isempty(pulses)
            if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
                listOfChannel = zeros(size(pulses.timestampsOn));
                for ii = 1:length(pulses.timestampsOn)
                    listOfChannel(ii) = ~isempty(pulses.timestampsOn{ii});
                end
                listOfChannel = find(listOfChannel);
            else
                listOfChannel = digitalChannelsList;
            end
        
            % CSD
            for mm = 1:length(listOfChannel)
                fprintf('Stimulus %3.i of %3.i \n',mm, length(listOfChannel)); %\n
                st = pulses.timestampsOn{listOfChannel(mm)};
                shanks = session.extracellular.electrodeGroups.channels;            
                shanks(excludeShanks) = [];
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:length(shanks)
                    lfp = getLFP(shanks{jj});
                    twin = 0.02;
                    [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DigitalCh #',num2str(listOfChannel(mm))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                saveas(gcf,['SummaryFigures\digitalPulsesCSD_ch',num2str(listOfChannel(mm)), '.png']);
%                 clear pulses listOfChannel
            end
        end

    catch
        warning('Error on CSD from digital inputs! ');
    end
end

% HIPPOCAMPAL LAYERS
if any(ismember(listOfAnalysis,{'hippocampalLayers'}))
    try
        [hippocampalLayers] = getHippocampalLayers('force',true,'promt',false);
    catch
        warning('Not possible to run getHippocampalLayers...')
    end
end
% PSTH RESPONSES ON DIGITAL AND ANALOG CH INPUTS
if any(ismember(listOfAnalysis,{'digitalPulses', 'analogPulses'}))
    try
        disp('Analog-in and/or digital-in PSTH...');
        % getting analog pulses channel
        if ~isempty(analogChannelsList)
            analogPulses = getAnalogPulses;
            if isempty(analogPulses)
                listOfAnalogChannel = [];
                if ~isempty(analogChannelsList)
                    warning('Analog channel list is empty, but some analog channel were specified for analysis!!')
                end
            else
                if ischar(analogChannelsList) && strcmpi(analogChannelsList,'all')
                    analogPulses = getAnalogPulses;
                    listOfAnalogChannel = unique(analogPulses.analogChannel);
                    clear analogPulses
                else
                    listOfAnalogChannel = analogChannelsList;
                end
            end
            clear analogPulses
        else
            listOfAnalogChannel = [];
        end
        
        % getting digital pulses channel
        if ~isempty(digitalChannelsList)
            digPulses = getDigitalIn;
            if isempty(digPulses)
                listOfDigitalChannel = [];
                if ~isempty(digitalChannelsList)
                    warning('Digital channel list is empty, but some digital channel were specified for analysis!!')
                end
            else
                if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
                    listOfDigitalChannel = zeros(size(digPulses.timestampsOn));
                    for ii = 1:length(digPulses.timestampsOn)
                        listOfDigitalChannel(ii) = ~isempty(digPulses.timestampsOn{ii});
                    end
                    listOfDigitalChannel = find(listOfDigitalChannel);
                    clear digPulses
                else
                    listOfDigitalChannel = digitalChannelsList;
                end
            end
            clear digPulses
        else
            listOfDigitalChannel = [];
        end
        
        optogeneticResponses = getOptogeneticResponse('analogChannelsList',listOfAnalogChannel,'digitalChannelsList', listOfDigitalChannel,'numRep',0);

    catch
        warning('Error on PSTH from digital and/or analog inputs! ');
    end
end

% DOWN-STATES
if any(ismember(listOfAnalysis,'downStates'))
    try
        disp('Slow-waves CSD and PSTH...');

        UDStates = detectUD;
        % CSD
        shanks = session.extracellular.electrodeGroups.channels;            
        shanks(excludeShanks) = [];
        twin = 0.2;
        evs = UDStates.timestamps.DOWN(:,1);
        figure
        set(gcf,'Position',[100 100 1400 600])
        for jj = 1:size(shanks,2)
            lfp = getLFP(shanks{jj},'noPrompts', true);
            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
            taxis = linspace(-twin,twin,size(csd.data,1));
            cmax = max(max(csd.data)); 
            subplot(1,size(shanks,2),jj);
            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DOWN-UP, Shank #',num2str(jj)),'FontWeight','normal'); 
            colormap jet; caxis([-cmax cmax]);
            hold on
            for kk = 1:size(lfpAvg.data,2)
                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
            end
        end
        saveas(gcf,'SummaryFigures\downUpStatesCSD.png');
        
        % PSTH
        psthUD = spikesPsth([],'eventType','slowOscillations','numRep',100,'min_pulsesNumber',10);
    catch
        warning('Error on Psth and CSD from down-states!');
    end
end

% RIPPLES
if any(ismember(listOfAnalysis,'ripples'))
    try
        disp('Ripples CSD and PSTH...');
        
        ripples = rippleMasterDetector;
        % CSD
        shanks = session.extracellular.electrodeGroups.channels;            
        shanks(excludeShanks) = [];
        twin = 0.1;
        evs = ripples.peaks;
        figure
        set(gcf,'Position',[100 100 1400 600])
        for jj = 1:size(shanks,2)
            lfp = getLFP(shanks{jj},'noPrompts', true);
            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
            taxis = linspace(-twin,twin,size(csd.data,1));
            cmax = max(max(csd.data)); 
            subplot(1,size(shanks,2),jj);
            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('Ripples, Shank #',num2str(jj)),'FontWeight','normal'); 
            colormap jet; caxis([-cmax cmax]);
            hold on
            for kk = 1:size(lfpAvg.data,2)
                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
            end
        end
        saveas(gcf,'SummaryFigures\ripplesCSD.png');
        
        % PSTH
        psthRipples = spikesPsth([],'eventType','ripples','numRep',100,'min_pulsesNumber',10);
        
    catch
        warning('Error on Psth and CSD from ripples! ');
    end
end

% THETA AND GAMMA PHASE MODULATION
if any(ismember(listOfAnalysis,'thetaModulation'))
    try
        thetaEpochs = detectThetaEpochs;
        computePhaseModulation();
    catch
        warning('It has not been possible to run theta and gamma mod code...');
    end
end

% TMAZEBEHAVIOUR AND LINEARMAZEBEHAVIOUR
if any(ismember(listOfAnalysis,'tMazeBehaviour')) || any(ismember(listOfAnalysis,'linearMazeBehaviour'))
   try 
        getSessionTracking('convFact',tracking_pixel_cm,'roiTracking','manual'); 
        if any(ismember(listOfAnalysis,'tMazeBehaviour'))
            getSessionArmChoice('task','alternation');
        end
        behaviour = getSessionLinearize;

        % PLACE CELLS SUMMARY
        spikes = loadSpikes('getWaveformsFromDat',false);
        firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
        placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps);
        firingTrialsMap = firingMapPerTrial;
   catch
       warning('It has not been possible to run the behaviour code...');
   end
end

cd(prevPath);
end